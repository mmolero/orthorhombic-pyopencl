"""  Class: ANI2D
     
     Miguel Molero (M. Molero)  & Ursula Iturraran-Viveros, 2012

	 miguel.molero@gmail,com
	 ursula.iturraran@gmail.com
	

"""

from   ANI2D_Classes  import *
import numpy		  as     np
import pyopencl		  as     cl




class ANI2D:
	
	def __init__(self, Image, Materials, Source, Transducer, Signal, SimModel):
	
		self.Im					  = np.float32(np.copy(SimModel.Im))
		self.MRI, self.NRI		  = np.shape(self.Im)
		self.Theta				  = Source.Theta
		self.Frequency            = Signal.Frequency
		self.StopSignal           = np.around( (2.0/self.Frequency)*(1/SimModel.dt) )
		self.Device			      = SimModel.Device
		
		
		self.InitCL(self.Device)
		self.MaterialSetup(Materials) 
		self.Init_Fields(Signal, SimModel)
		self.ReceiverSetup(Image,Source,Transducer, SimModel)
		self.StaggeredProp()
		self.applyABC(Materials, SimModel)
		self.n		   = 0
		
		self.Init_Fields_CL(SimModel)
		
		
		
	def InitCL(self, DEVICE="GPU"):
		
		try:
			for platform in cl.get_platforms():
				for device in platform.get_devices():
					if cl.device_type.to_string(device.type)== DEVICE:
						my_device =	 device
						print my_device.name, "	 ", cl.device_type.to_string(my_device.type)
						
		except:
			my_device = cl.get_platforms()[0].get_devices()
			print my_device.name, "	 ", cl.device_type.to_string(my_device.type)
			
			
						
		self.ctx   = cl.Context([my_device])
		self.queue = cl.CommandQueue(self.ctx)
		self.mf	   = cl.mem_flags
		
		
	def MaterialSetup(self, Materials):

		NumeroMat  = Materials.NumMat
		#Vacuum Condition if some medium is air 
		for n in range(0,NumeroMat):
			if	Materials.rho[n] < 2.0:
				Materials.rho[n] = 10e23
				Materials.c11[n] = 1e-20
				Materials.c12[n] = 1e-20
				Materials.c22[n] = 1e-20
				Materials.c44[n] = 1e-20

		self.rho	= np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.c11	= np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.c12	= np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.c22	= np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.c44	= np.ones((self.MRI,self.NRI) ,dtype=np.float32)*1e-30

		for i in range(0,self.MRI):
			for n in range(0,NumeroMat):
				ind =  np.nonzero(self.Im[i,:] == Materials.gray[n])
				self.rho[i,ind] = Materials.rho[n]
				self.c11[i,ind] = Materials.c11[n]
				self.c12[i,ind] = Materials.c12[n]
				self.c22[i,ind] = Materials.c22[n]
				self.c44[i,ind] = Materials.c44[n]
				if self.c44[i,ind].any() == 0:
					self.c44[i,ind] = 1e-30
					

		for i in range(0,self.MRI):
			ind = np.nonzero(self.Im[i,:] == 255.0)
			self.rho[i,ind]	 = Materials.rho[0]
			self.c11[i,ind]	 = Materials.c11[0]
			self.c12[i,ind]	 = Materials.c12[0]
			self.c22[i,ind]	 = Materials.c22[0]
			self.c44[i,ind]	 = Materials.c44[0]
			if self.c44[i,ind].any() == 0:
				self.c44[i,ind] = 1e-30
			
	
	def ConfigFuente(self, Image, Source, Transducer, SimModel):

		if len(Transducer) == 2: 
			
			self.MoreTransducer = True
			self.insp     = self.setInspection(Image,Source,Transducer[0],SimModel)
			
			self.XL       = np.copy(self.insp.XL)
			self.YL       = np.copy(self.insp.YL)
			self.IR       = np.copy(self.insp.IR)
			self.insp_R   = self.setInspection(Image,Source,Transducer[1],SimModel)
			
		else:
			self.MoreTransducer = False
			self.insp     = self.setInspection(Image,Source, Transducer, SimModel)
			self.XL       = np.copy(self.insp.XL)
			self.YL       = np.copy(self.insp.YL)
			self.IR       = np.copy(self.insp.IR)
			
			
			
	def setInspection(self,Image,Source, Transducer, SimModel):
	
		insp = Inspection()
		
		D_T = (self.MRI-2.)/2.
		x2	= self.MRI/2.  + (D_T - SimModel.TapG)*np.sin(Source.Theta)
		y2	= self.NRI/2.  + (D_T - SimModel.TapG)*np.cos(Source.Theta)

		X0	= self.MRI/2.
		Y0	= self.NRI/2.

		Transducer.SizePixel =  np.around( 0.5 * Image.Pixel_mm * Transducer.Size	 * float(self.NRI) / Image.N_pml )
		
		insp.setTransmisor(Source,Transducer,x2,y2,X0,Y0)
		insp.addOffset(Image, Transducer, self.NRI)
		insp.addBorderOffset(Image, Transducer, self.MRI)
		
		return insp
		
		
	def StaggeredProp(self):

		BXtemp		 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		BYtemp		 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.BX		 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.BY		 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.C11	 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.C12	 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.C22	 = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.C44	 = np.zeros((self.MRI,self.NRI),dtype=np.float32)

		BXtemp[:,:]	 =	1.0/self.rho[:,:]
		BYtemp[:,:]	 =	1.0/self.rho[:,:]

		self.C11	 =	 np.copy(self.c11)
		self.C12	 =	 np.copy(self.c12)
		self.C22	 =	 np.copy(self.c22)

		self.BX[:-2,:]	   =  0.5*( BXtemp[1:-1,:] + BXtemp[:-2,:] )
		self.BX[ -2,:]	   =  np.copy(BXtemp[-2,:])

		self.BY[:,:-2]	   =  0.5*( BYtemp[:,1:-1]	 + BYtemp[:,:-2]  )
		self.BY[:, -2]	   =  np.copy(BYtemp[:,-2])

		self.C44[:-2,:-2]  = 4./(  (1./self.c44[:-2,:-2] ) +  (1./self.c44[1:-1,:-2]) +	 (1./self.c44[:-2,1:-1] ) +	 (1./self.c44[1:-1,1:-1] )	)

			

	def ConfigAirBoundary(self):

		indx,indy = np.nonzero(self.Im == 255)
		self.BX[indx,indy]	 = 0.0
		self.BY[indx,indy]	 = 0.0
		self.C11[indx,indy]	 = 0.0
		self.C12[indx,indy]	 = 0.0
		self.C22[indx,indy]	 = 0.0
		self.C44[indx,indy]	 = 0.0

		
		
	def Init_Source(self, Signal,t):
		self.input_source  = Signal.generate(t)

	def Init_Fields(self, Signal, SimModel):
		
		self.input_source  = Signal.generate(SimModel.t)
		
		self.vx			   = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.vy			   = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.Txx		   = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.Txy		   = np.zeros((self.MRI,self.NRI),dtype=np.float32)
		self.Tyy		   = np.zeros((self.MRI,self.NRI),dtype=np.float32)
			
		self.SV			   = np.zeros( (self.MRI,self.NRI), dtype=np.float32)
		self.ABC		   = np.zeros( (self.MRI,self.NRI), dtype=np.float32)
		
		
	
	def Init_Fields_CL(self, SimModel):
		
		self.Txx_buf	   = cl.Buffer(self.ctx, self.mf.READ_WRITE	| self.mf.COPY_HOST_PTR, hostbuf=self.Txx)
		self.Tyy_buf	   = cl.Buffer(self.ctx, self.mf.READ_WRITE	| self.mf.COPY_HOST_PTR, hostbuf=self.Tyy)
		self.Txy_buf	   = cl.Buffer(self.ctx, self.mf.READ_WRITE	| self.mf.COPY_HOST_PTR, hostbuf=self.Txy)
		self.vx_buf		   = cl.Buffer(self.ctx, self.mf.READ_WRITE	| self.mf.COPY_HOST_PTR, hostbuf=self.vx)
		self.vy_buf		   = cl.Buffer(self.ctx, self.mf.READ_WRITE	| self.mf.COPY_HOST_PTR, hostbuf=self.vy)	
		
		
		self.BX_1		   = np.copy(self.BX)
		self.BY_1		   = np.copy(self.BY)
		self.C11_1		   = np.copy(self.C11)
		self.C12_1		   = np.copy(self.C12)
		self.C22_1		   = np.copy(self.C22)
		self.C44_1		   = np.copy(self.C44)
		
		self.ConfigSource()

		self.ABC_buf	   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.ABC)
		self.BX_buf		   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.BX)
		self.BY_buf		   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.BY)
		self.C11_buf	   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.C11)
		self.C12_buf	   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.C12)
		self.C22_buf	   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.C22)
		self.C44_buf	   = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.C44)

		
		self.NX	           = np.size(self.XL,0)
		
		self.XLL	       = np.copy(np.int32(self.XL[:,0]))	 
		self.YLL	       = np.copy(np.int32(self.YL[:,0]))
		
		self.XL_buf	       = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.XLL)
		self.YL_buf	       = cl.Buffer(self.ctx, self.mf.READ_ONLY	| self.mf.COPY_HOST_PTR, hostbuf=self.YLL)
		
		
		self.dtx		   = np.float32(SimModel.dt/SimModel.dx)
		self.dtdxx		   = np.float32(SimModel.dt/(SimModel.dx**2))
		
	
			
		self.program	    = cl.Program( self.ctx, self.ANI2D_Kernel() ).build()
		
		
		
	def applyABC(self,Materials, SimModel): 
		
		APARA = 0.015										
		for i in range(0,self.MRI):													 
			for j in range(0,self.NRI):
				if	 i < SimModel.TapG:											 
					self.ABC[i,j] = np.exp(-((APARA*(SimModel.TapG-i))**2))				 
				elif j < SimModel.TapG:									  
					self.ABC[i,j] = np.exp(-((APARA*(SimModel.TapG-j))**2))			   
				elif i > (self.MRI-SimModel.TapG+1):							 
					self.ABC[i,j] = np.exp(-((APARA*(i-self.MRI+SimModel.TapG-1))**2))					
				elif j > (self.NRI-SimModel.TapG+1):							 
					self.ABC[i,j] = np.exp(-((APARA*(j-self.NRI+SimModel.TapG-1))**2))				   
				else:															
					self.ABC[i,j] = 1.0

		
	def ConfigSource(self):

		NX		= np.size(self.XL,0)
		IT		= -1
		IT1		= -2
		
		
		for ITT	in range(-2,3):
			for IT in range(-3,0):
				for m in range(0,NX):

					xl	= int(self.XL[m,0])
					yl	= int(self.YL[m,0])
   
					self.BX[xl+IT,yl+ITT]	= 0.0
					self.BY[xl+IT,yl+ITT]	= 0.0
					self.C11[xl+IT,yl+ITT]	= 0.0
					self.C12[xl+IT,yl+ITT]	= 0.0
					self.C22[xl+IT,yl+ITT]	= 0.0
					self.C44[xl+IT,yl+ITT]	= 0.0
					
			
		
	def ReceiverSetup(self, Image,Source,Transducer, SimModel):
	
		self.ConfigFuente(Image,Source,Transducer,SimModel)
		
		if self.MoreTransducer:
			self.receiver_signals_R = np.zeros(( SimModel.Ntime, np.size(self.IR,1)) ,dtype=np.float32)
			
		self.receiver_signals = np.zeros(( SimModel.Ntime, np.size(self.IR,1) ) ,dtype=np.float32)
		
	
	def AgainReceiverSetup(self, Image,Source,Transducer, SimModel,Ntime):

		self.ConfigFuente(Image,Source,Transducer,SimModel)

		if self.MoreTransducer:
			self.receiver_signals_R = np.zeros(( Ntime, np.size(self.IR,1)) ,dtype=np.float32)

		self.receiver_signals = np.zeros(( Ntime, np.size(self.IR,1) ) ,dtype=np.float32)
	
	
	def receivers(self):
	
		cl.enqueue_copy(self.queue, self.Txx, self.Txx_buf)
		if self.MoreTransducer:
			self.receiver_signals_R[self.n,:] = self.insp_R.SetReception(self.Txx)
			
		self.receiver_signals[self.n,:] = self.insp.SetReception(self.Txx)

		
			
	def RunSim(self):
		
		self.Run_Global()
			
		#air backing for transducer emitter exits when there be signal
		if self.n == self.StopSignal:
			print "iter stop signal", self.StopSignal
			self.BX		   = np.copy(self.BX_1)
			self.BY		   = np.copy(self.BY_1)
			self.C11	   = np.copy(self.C11_1)
			self.C12	   = np.copy(self.C12_1)
			self.C22	   = np.copy(self.C22_1)
			self.C44	   = np.copy(self.C44_1)
			cl.enqueue_copy(self.queue, self.BX_buf, self.BX)
			cl.enqueue_copy(self.queue, self.BY_buf, self.BY)
			cl.enqueue_copy(self.queue, self.C11_buf, self.C11)
			cl.enqueue_copy(self.queue, self.C12_buf, self.C12)
			cl.enqueue_copy(self.queue, self.C22_buf, self.C22)
			cl.enqueue_copy(self.queue, self.C44_buf, self.C44)
				
		self.receivers()
		
	
						
		
	def Run_Global(self):
		
										
		self.program.VelX_ANI2D(self.queue, (self.NRI,self.MRI,), None,
								    self.Txx_buf,
									self.Txy_buf,
									self.vx_buf,
									self.BX_buf,
									self.ABC_buf)
											
		self.program.VelY_ANI2D(self.queue, (self.NRI,self.MRI,), None,
									self.Txy_buf,
									self.Tyy_buf,
									self.vy_buf,
									self.BY_buf,
									self.ABC_buf)
											 
										
		self.program.StressX_ANI2D(self.queue, (self.NRI,self.MRI,), None,
									self.Txx_buf,
								    self.vx_buf,
									self.vy_buf,
									self.C11_buf,
									self.C12_buf,
									self.ABC_buf)
											
		self.program.StressY_ANI2D(self.queue, (self.NRI,self.MRI,), None,
							        self.Tyy_buf,
									self.vx_buf,
									self.vy_buf,
									self.C11_buf,
									self.C12_buf,
									self.C22_buf,
								    self.ABC_buf)
								
		self.program.StressXY_ANI2D(self.queue, (self.NRI,self.MRI,), None,
									self.Txy_buf,
									self.vx_buf,
									self.vy_buf,
									self.C44_buf,
									self.ABC_buf)
										

		
		y  = np.float32(self.input_source[self.n])
		
		self.program.Source_ANI2D( self.queue, (self.NX,), None,
		                                   	 self.Txx_buf,
											 self.XL_buf,
											 self.YL_buf, 
											 y)
											
	
	def RunGL(self, step=50):

		if self.n % step==0:
			cl.enqueue_copy(self.queue, self.vx, self.vx_buf)
			cl.enqueue_copy(self.queue, self.vy, self.vy_buf)
			self.SV	 = np.sqrt(self.vx**2 + self.vy**2 )
			self.SV	 = 20.*np.log10((np.abs(self.SV)/np.max(np.abs(self.SV+1e-40))) + 1e-40)



	def ANI2D_Kernel(self):

		macro	 = """#define MRI		%s
					  #define NRI		%s
					  #define ind(i, j)	 ( ( (i)*NRI) + (j) )
					  #define dtx		%gf
					  #define dtdxx		%gf
					 

		"""%( str(self.MRI), str(self.NRI), self.dtx, self.dtdxx)
		kernel_source = macro + """



           __kernel void Source_ANI2D(	__global float *Txx, 
                                        __global const int *XL, __global const int *YL,  
                                          const float source){

			     uint m =  get_global_id(0);
			     
		         int ix = XL[m];
			     int iy = YL[m];
			     Txx[ind(ix,iy)] -= (source*dtdxx);

           }

		   __kernel void VelX_ANI2D( __global float *Txx,	__global float *Txy,	   
								     __global float *vx,	__global const float *BX,	
									 __global const float *ABC){
									
			  

				   float    W80 = 1.19628906E+00f;
				   float    W81 = 7.97526017E-02f;
				   float    W82 = 9.57031269E-03f;
				   float    W83 = 6.97544659E-04f;

				   int j = get_global_id(0);
				   int i = get_global_id(1);


					if(  i > 3 &&  i <MRI-4  && j >3 && j<NRI-4 ){

					vx[ind(i,j)]  +=  ( (  dtx  * BX[ind(i,j)] ) *( 
					                    (  Txy[ind(i,j)]   -Txy[ind(i,j-1)]   )*W80
									  - (  Txy[ind(i,j+1)] -Txy[ind(i,j-2)]   )*W81
									  + (  Txy[ind(i,j+2)] -Txy[ind(i,j-3)]   )*W82
									  - (  Txy[ind(i,j+3)] -Txy[ind(i,j-4)]   )*W83    +  
					
					                    (  Txx[ind(i+1,j)] -Txx[ind(i,j)]     )*W80
									  - (  Txx[ind(i+2,j)]- Txx[ind(i-1,j)]   )*W81
									  + (  Txx[ind(i+3,j)]- Txx[ind(i-2,j)]   )*W82
									  - (  Txx[ind(i+4,j)]- Txx[ind(i-3,j)]   )*W83   ) );
					
					}

			        barrier(CLK_GLOBAL_MEM_FENCE);
              
			        vx[ind(i,j)]	  *= ABC[ind(i,j)];
			         
			  
		
			 }

			__kernel void VelY_ANI2D( __global float *Txy,	__global float *Tyy,
									  __global float *vy, 	__global const float *BY,
								      __global const float *ABC){


					   float    W80 = 1.19628906E+00f;
					   float    W81 = 7.97526017E-02f;
					   float    W82 = 9.57031269E-03f;
					   float    W83 = 6.97544659E-04f;

					   int j = get_global_id(0);
					   int i = get_global_id(1);


					   if(  i > 3 &&  i <MRI-4  && j >3 && j<NRI-4 ){

						vy[ind(i,j)]  +=  ( (  dtx * BY[ind(i,j)] ) *(   
						                    (  Txy[ind(i,j)]   -Txy[ind(i-1,j)]   )*W80
										  - (  Txy[ind(i+1,j)] -Txy[ind(i-2,j)]   )*W81
										  + (  Txy[ind(i+2,j)] -Txy[ind(i-3,j)]   )*W82
										  - (  Txy[ind(i+3,j)] -Txy[ind(i-4,j)]   )*W83   +

										    (  Tyy[ind(i,j+1)] -Tyy[ind(i,j)]     )*W80
										  - (  Tyy[ind(i,j+2)]- Tyy[ind(i,j-1)]   )*W81
										  + (  Tyy[ind(i,j+3)]- Tyy[ind(i,j-2)]   )*W82
										  - (  Tyy[ind(i,j+4)]- Tyy[ind(i,j-3)]   )*W83   )  );

	                    }


				        barrier(CLK_GLOBAL_MEM_FENCE);
				
				        vy[ind(i,j)]	  *= ABC[ind(i,j)];	  


				 }




		  __kernel void StressX_ANI2D( __global float *Txx,	     
								       __global float *vx,		  __global float *vy,
								       __global const float *C11,  __global const float *C12, 
								       __global const float *ABC) {
				
				
				
			
		     float    W80 = 1.19628906E+00f;
			 float    W81 = 7.97526017E-02f;
			 float    W82 = 9.57031269E-03f;
			 float    W83 = 6.97544659E-04f;

			 int j = get_global_id(0);
			 int i = get_global_id(1);
			
	              
			 if(  i > 3 &&  i <MRI-4  && j >3 && j<NRI-4 ){
			
			 
	                 Txx[ind(i,j)] += ( (dtx * C11[ind(i,j)] ) *  ( (  vx[ind(i,j)]   -vx[ind(i-1,j)] )*W80
											                      - (  vx[ind(i+1,j)] -vx[ind(i-2,j)] )*W81
											                      + (  vx[ind(i+2,j)] -vx[ind(i-3,j)] )*W82
											                      - (  vx[ind(i+3,j)] -vx[ind(i-4,j)] )*W83 ) +
	
	                                    (dtx * C12[ind(i,j)] ) *  ( (  vy[ind(i,j)]   -vy[ind(i,j-1)] )*W80
																  - (  vy[ind(i,j+1)] -vy[ind(i,j-2)] )*W81
																  + (  vy[ind(i,j+2)] -vy[ind(i,j-3)] )*W82
																  - (  vy[ind(i,j+3)] -vy[ind(i,j-4)] )*W83 ) );

            }          

									
			 barrier(CLK_GLOBAL_MEM_FENCE);
           
			 Txx[ind(i,j)]	*= ABC[ind(i,j)];
			 
		      
		 
		  }
							
							
		 __kernel void StressY_ANI2D(  __global float *Tyy,
								      __global float *vx,		  __global float *vy,
								      __global const float *C11,  __global const float *C12, 
								      __global const float *C22,  
								      __global const float *ABC) {




		     float    W80 = 1.19628906E+00f;
			 float    W81 = 7.97526017E-02f;
			 float    W82 = 9.57031269E-03f;
			 float    W83 = 6.97544659E-04f;

			 int j = get_global_id(0);
			 int i = get_global_id(1);


			 if(  i > 3 &&  i <MRI-4  && j >3 && j<NRI-4 ){

				     Tyy[ind(i,j)] += ( (dtx * C12[ind(i,j)] ) *  ( (  vx[ind(i,j)]   -vx[ind(i-1,j)] )*W80
																  - (  vx[ind(i+1,j)] -vx[ind(i-2,j)] )*W81
																  + (  vx[ind(i+2,j)] -vx[ind(i-3,j)] )*W82
																  - (  vx[ind(i+3,j)] -vx[ind(i-4,j)] )*W83 ) +

									    (dtx * C22[ind(i,j)] ) *  ( (  vy[ind(i,j)]   -vy[ind(i,j-1)] )*W80
																  - (  vy[ind(i,j+1)] -vy[ind(i,j-2)] )*W81
																  + (  vy[ind(i,j+2)] -vy[ind(i,j-3)] )*W82
																  - (  vy[ind(i,j+3)] -vy[ind(i,j-4)] )*W83 ) );

            }          


			 barrier(CLK_GLOBAL_MEM_FENCE);
			
			 Tyy[ind(i,j)]	*= ABC[ind(i,j)];
			 


		  }
				
				
		 __kernel void StressXY_ANI2D( __global float *Txy,
								       __global float *vx,		  __global float *vy,
								       __global const float *C44,
								       __global const float *ABC) {


		     float    R80 = 1.19628906E+00f;
			 float    R81 = 7.97526017E-02f;
			 float    R82 = 9.57031269E-03f;
			 float    R83 = 6.97544659E-04f;

			 int j = get_global_id(0);
			 int i = get_global_id(1);

			if(  i > 3 &&  i <MRI-4  && j >3 && j<NRI-4 ){		

				     Txy[ind(i,j)] +=  (  (dtx  * C44[ind(i,j)] )*( (  vx[ind(i,j+1)] -vx[ind(i,j)]   )*R80
																  - (  vx[ind(i,j+2)]- vx[ind(i,j-1)] )*R81
																  + (  vx[ind(i,j+3)]- vx[ind(i,j-2)] )*R82
																  - (  vx[ind(i,j+4)]- vx[ind(i,j-3)] )*R83  +

																    (  vy[ind(i+1,j)] -vy[ind(i,j)]   )*R80
																  - (  vy[ind(i+2,j)]- vy[ind(i-1,j)] )*R81
																  + (  vy[ind(i+3,j)]- vy[ind(i-2,j)] )*R82
																  - (  vy[ind(i+4,j)]- vy[ind(i-3,j)] )*R83  ) );   
			 }



			 barrier(CLK_GLOBAL_MEM_FENCE);

			 Txy[ind(i,j)]	*= ABC[ind(i,j)];


		  }
							
													
								
		"""
		return kernel_source
		
		



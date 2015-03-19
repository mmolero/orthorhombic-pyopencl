"""  Library of used classes

     Miguel Molero (M. Molero)  & Ursula Iturraran-Viveros, 2012

	 miguel.molero@gmail,com
	 ursula.iturraran@gmail.com
	


"""


import   numpy  as  np
from	 math				import sin, cos, sqrt, pi
from	 scipy.misc			import imresize
from	 scipy.misc.pilutil import imrotate
from	 scipy.weave		import inline
from	 scipy.io			import savemat
#from	 PIL                import Image



class NewImage:

	def __init__(self, Width=40, Height=40,Pixel_mm=10,GrayScale=0,SPML=False):
		self.Width		 = Width
		self.Height		 = Height
		self.Pixel_mm	 = Pixel_mm
		self.GrayScale	 = GrayScale
			
		self.M			 = int(self.Height * self.Pixel_mm)
		self.N			 = int(self.Width  * self.Pixel_mm)
		self.I			 = np.ones((self.M,self.N),dtype=np.uint8)*GrayScale
		
	def createABC(self,Tap):

		self.Tap		 = Tap
	
		TP				 = round(Tap* self.Pixel_mm )
		self.M_pml		 = int( self.M	 + 2*TP )
		self.N_pml		 = int( self.N	 + 2*TP )

		self.Itemp		 = 255.0*np.ones((self.M_pml,self.N_pml),dtype=np.uint8)
		self.Itemp[TP : self.M_pml-TP, TP : self.N_pml-TP] = np.copy(self.I)

#####################################################################################
		
	
def addRectangle(NewImage, W_0, H_0, W, H, GrayScale, Theta):

	a  = int(H*NewImage.Pixel_mm/2) 
	b  = int(W*NewImage.Pixel_mm/2) 
	for	 x in  range(-a,a+1):
		for y in range(-b,b+1):
			tempX = round (x + H_0*NewImage.Pixel_mm)
			tempY = round (y + W_0*NewImage.Pixel_mm)
			NewImage.I[tempX,tempY] = GrayScale

	if Theta != 0:
		NewImage.I = imrotate(NewImage.I,Theta,interp='nearest')

	return NewImage.I


###################################################################################

class Material:

	def __init__(self, name="Water",rho=1000,c11=2.19e9,c12=0.0,c22=0.0,c44=0.0,gray=0):

		self.name	=  name
		self.rho	=  rho
		self.c11	=  c11
		self.c12	=  c12
		self.c22	=  c22
		self.c44	=  c44
		self.VL		=  sqrt( c11/rho )
		self.VT		=  sqrt( c44/rho )
		self.gray	=  gray


class MaterialSet:

	def __init__(self):

		self.rho	= []
		self.c11	= []
		self.c12	= []
		self.c22	= []
		self.c44	= []
		self.name	= []
		self.gray	= []
		self.NumMat = 0
		self.VL		= []
		self.VT		= []

	def append(self, material):

		self.name.append(material.name)
		self.rho.append(material.rho)
		self.c11.append(material.c11)
		self.c12.append(material.c12)
		self.c22.append(material.c22)
		self.c44.append(material.c44)
		self.gray.append(material.gray)
		self.VL.append(material.VL)
		self.VT.append(material.VT)	 

		self.NumMat += 1

	def remove(self,material):

		self.name.remove(material.name)
		self.rho.remove(material.rho)
		self.c11.remove(material.c11)
		self.c12.remove(material.c12)
		self.c22.remove(material.c22)
		self.c44.remove(material.c44)
		self.gray.remove(material.gray)
		self.VL.remove(material.VL)
		self.VT.remove(material.VT)	 

		self.NumMat -= 1
		if self.NumMat < 0:
			self.NumMat=0

	def delete(self):

		self.name.pop()
		self.rho.pop()
		self.c11.pop()
		self.c12.pop()
		self.c22.pop()
		self.c44.pop()
		self.gray.pop()
		self.VL.pop()
		self.VT.pop()	 

		self.NumMat -= 1
		if self.NumMat < 0:
			self.NumMat=0


#############################################################################################


class Source:

	def __init__(self,LaunchType = 'Dual'):

		self.LaunchType		= LaunchType
		self.Theta			= 0

		if	 self.LaunchType == 'PulseEcho':
			self.pulseEcho()

		elif self.LaunchType == 'Transmission':
			self.transmission()
			
		elif self.LaunchType == "Dual":
			self.dual()


	def pulseEcho(self):
		self.Theta = [270*pi/180.0, 270.0*pi/180.0]
		#print self.Theta

	def transmission(self):
		self.Theta = [270.0*pi/180.0, 90.0*pi/180.0]
		#print self.Theta

	def dual(self):
		self.Theta = [270.0*pi/180.0, 270*pi/180.0, 90.0*pi/180.0]


################################################################################################

class Transducer:

	def __init__(self, Size = 10, Offset=0, BorderOffset=0, Location=0):
		# Location = 0 => Top
		self.Size		  = Size
		self.Offset		  = Offset
		self.BorderOffset = BorderOffset
		self.SizePixel	  = []
		self.Location	  = Location
		

################################################################################################

class Signal:

	def __init__(self, Option=0, Amplitude=1, Frequency=1, N_Cycles=1):

		self.Amplitude = Amplitude
		self.Frequency = Frequency
		self.Option	   = Option
	
		if Option == 0:
			self.name  = "Raised Cosine Pulse"
		
	def generate(self,t):

		if self.Option == 0:
			return RaisedCosinePulse(t, self.Frequency, self.Amplitude)
	
	def saveSignal(self,t):	

		self.time_signal  = self.generate(t)


################################################################################################


def RaisedCosinePulse(t, Freq, Amplitude):

	N  = np.size(t,0)
	P  = np.zeros((N,),dtype=np.float32)
	for m in range(0,N):
		if t[m] <= 2.0/Freq :
			P[m] = Amplitude *(1-cos(pi*Freq*t[m]))*cos(2*pi*Freq*t[m])

	return P



################################################################################################

class Inspection:

	def __init__(self):

		self.Theta	= []
		self.XL		= []
		self.YL		= []
		self.IR		= []


	def setTransmisor(self, source, transducer, x2, y2, X0, Y0):

		self.Theta	= source.Theta

		Ntheta		= np.size(self.Theta,0)
		NXL			= int(2*transducer.SizePixel)

		xL			= np.zeros((NXL,),dtype=np.float32)
		yL			= np.zeros((NXL,),dtype=np.float32)

		for m in range(0,Ntheta):

			if np.abs(np.cos(self.Theta[m])) < 1e-5:
				yL = np.linspace(y2[m]-transducer.SizePixel,y2[m]+transducer.SizePixel,num=NXL, endpoint=True)
				xL[:] = x2[m]*np.ones((NXL,),dtype=np.float32)


			elif np.abs(np.cos(self.Theta[m])) == 1:
				xL[:] = np.linspace(x2[m]-transducer.SizePixel, x2[m]+transducer.SizePixel,num=NXL, endpoint=True)
				yL[:] = y2[m] - ( (x2[m]-X0 )/( y2[m]-Y0 ) )*( xL[:]-x2[m] )

			else:
				xL[:] = np.linspace(x2[m]-(transducer.SizePixel*np.abs(np.cos(self.Theta[m]))),x2[m]+(transducer.SizePixel*np.abs(np.cos(self.Theta[m]))), num=NXL, endpoint=True )
				yL[:] = y2[m] - ( (x2[m]-X0 )/( y2[m]-Y0 )	)*( xL[:]-x2[m] )

			if m==0:
				self.XL		= np.zeros((np.size(xL,0),Ntheta),dtype=np.float32)
				self.YL		= np.zeros((np.size(xL,0),Ntheta),dtype=np.float32)


			self.XL[:,m]  = (np.around(xL[:]))
			self.YL[:,m]  = (np.around(yL[:]))



	def addOffset(self, image, transducer, NRI):

		NXL	   = np.size(self.XL,0)
		Ntheta = np.size(self.Theta,0)

		self.YL +=	 (np.around(transducer.Offset * image.Pixel_mm * NRI / float(image.N_pml)))

		self.IR		 = np.zeros((Ntheta,Ntheta-1),dtype=np.float32)
		B			 = range(0,Ntheta)
		self.IR[:,0] = np.int32(B[:])

		for i in range(1,Ntheta-1):
			B  = np.roll(B,-1)
			self.IR[:,i] = np.int32(B)

	def addBorderOffset(self, image, transducer, MRI):

		self.XL[:,0] += (np.around(transducer.BorderOffset * image.Pixel_mm * MRI / float(image.M_pml)))
		self.XL[:,1] -= (np.around(transducer.BorderOffset * image.Pixel_mm * MRI / float(image.M_pml)))


	def SetReception(self,T):

		ReceptorX = (self.XL)
		ReceptorY = (self.YL)
		M,N		  = np.shape(ReceptorX)
		temp  = np.zeros((M,N-1),dtype=np.float32)

		for	 mm	 in range(0,M):
			for ir in  range(0,N-1):

				temp[mm,ir]	  =	 T[ int(ReceptorX[ mm,int(self.IR[ir+1,0]) ] ) , int(ReceptorY[ mm,int(self.IR[ir+1,0]) ]) ]

		return np.mean(temp,0)


###########################################################################################################

class SimulationModel:
	
	def __init__(self,TimeScale=1, MaxFreq=2, PointCycle=10, SimTime=50):

		self.TimeScale	= TimeScale
		self.MaxFreq	= MaxFreq	  # MHz
		self.PointCycle = PointCycle
		self.SimTime	= SimTime	  # microseconds
		
		
	def jobParameters(self,materials):
		
		indVL = np.array(materials.VL) > 400
		indVT = np.array(materials.VT) > 400
		VL	  = np.array(materials.VL)[indVL]
		VT	  = np.array(materials.VT)[indVT]
		V	  = np.hstack( (VL, VT) )
		
		self.dx = np.float32( np.min([V]) / (self.PointCycle*self.MaxFreq) )
		
		
		self.dt = self.TimeScale * np.float32( 0.5* self.dx / (	 np.max([V]) ) )

		self.Ntime   		 = round(self.SimTime/self.dt)
		self.t				 = self.dt*np.arange(0,self.Ntime)
	
	
	def createNumericalModel(self,image):  
		
		#Spatial Scale
		Mp			      =	 image.M_pml*1e-3/image.Pixel_mm/self.dx
		self.Rgrid	      =	 Mp/image.M_pml
		
		self.TapG	      =	 np.around(image.Tap * self.Rgrid * image.Pixel_mm)
		self.Im		      =	 imresize(image.Itemp, self.Rgrid, interp='nearest')
		self.MRI,self.NRI =	 np.shape(self.Im)
	
		
	def initReceivers(self):
		
		self.receiver_signals   = []
		self.receiver_signals_R = []
		
	
	def setDevice(self,Device):
	
		self.Device = Device
		
			
		
	def rotatedModel(self, image, theta):
		Theta = theta*(180.0/pi) - 270.
		
		I            = imrotate(image.I,Theta,interp='nearest')
		#if np.shape(I) != np.shape(image.I):
		#	I = I.transpose()
		
		TP			 = int( image.Tap* image.Pixel_mm )
		M_pml		 = int( image.M	 + 2*TP )
		N_pml		 = int( image.N	 + 2*TP )
		Itemp		 = 255.0*np.ones((M_pml,N_pml),dtype=np.uint8)
		Itemp[TP : M_pml-TP, TP : N_pml-TP] = np.copy(I)
	
	
		Mp			      =	 M_pml*1e-3/image.Pixel_mm/self.dx
		self.Rgrid	      =	 Mp/M_pml

		self.TapG	      =	 np.around(image.Tap * self.Rgrid * image.Pixel_mm)
		self.Im		      =	 imresize(Itemp, self.Rgrid, interp='nearest')
		self.MRI,self.NRI =	 np.shape(self.Im)
			
			
		
	




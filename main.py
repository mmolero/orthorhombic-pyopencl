"""
Main program

Miguel Molero (M. Molero)  & Ursula Iturraran-Viveros, 2012

miguel.molero@gmail,com
ursula.iturraran@gmail.com

"""

from   ANI2D_Classes import *
from   ANI2D         import *

import sys
import time
from   scipy.io      import savemat


file_data      = 'test'  # it will be stored in matlab format

Width          = 75      # Dimension in mm
Height         = 100     # Dimension in mm
Tap            = 4       # Dimension in mm, Absorbing layer

Width_Sample   = 35      # Dimension in mm, Width Plate
Height_Sample  = 50      # Dimension in mm, Height Plate

Frequency      = 500e3   # operating frequency
SimTime        = 100e-6; # time simulation
device         = 'GPU'   # CPU or GPU

# rotation angles
THETA = np.linspace(-50,50,101, endpoint=True)*pi/180.0



if __name__ == '__main__':

	#import   matplotlib
	from  matplotlib.pylab   import *

	#Create Scenario
	Scenario    = NewImage(Width,Height,Pixel_mm = 10, GrayScale=0)
	Scenario.I  = addRectangle(Scenario, int(Width*0.5), int(Height*0.5), Width_Sample, Height_Sample, GrayScale=60, Theta=0)
	
	#Material Setup
	M1          = Material( name="M1",rho=1000,c11=2.19e9,c12=2.19e9,c22=2.19e9,c44=1e-30,gray=0)
	M2          = Material( name="M2",rho=2648,c11=15.9e9,c12=7.0e9,c22=15.9e9,c44=3.4e9,gray= 60)
	materials   = MaterialSet()
	materials.append(M1)
	materials.append(M2)
		
	#Absorbing conditions
	Scenario.createABC(Tap)
	
	#Create Inspection Type, Dual = Reflection & Transmission Inspection
	source       = Source(LaunchType = 'Dual')     
	emitter      = Transducer(Size=50,  Offset=0)
	receiver     = Transducer(Size=2.0, Offset=0)
	
	#Transducers
	transducers  = [emitter, receiver]
	
	#Input Signal
	waveform     = Signal(Amplitude=1, Frequency=Frequency)
	
	
	sim_model    = SimulationModel(TimeScale=1, MaxFreq=2*Frequency, PointCycle=10, SimTime=SimTime)
	sim_model.jobParameters(materials)
	sim_model.createNumericalModel(Scenario)
	sim_model.initReceivers()
	sim_model.setDevice(device)
	print sim_model.MRI, sim_model.NRI

	
	start = time.time()
	
	#Transmission and Reflection Matrices 
	MatrizT = np.zeros((sim_model.Ntime,np.size(THETA,0)),dtype=np.float32)
	MatrizR = np.zeros((sim_model.Ntime,np.size(THETA,0)),dtype=np.float32)
	
	nt = 0
	for theta in THETA:
		#rotation of scenario
		if theta != 0:
			sim_model.rotatedModel(Scenario, theta)
		#Create Simulation Object
		FD = ANI2D(Scenario, materials, source, transducers, waveform, sim_model)

		#loop computing wave propagation
		while (FD.n < sim_model.Ntime):
			FD.RunSim()
			if FD.n % 500==0:
				print FD.n, " of iterations: ", sim_model.Ntime, " rotations: ", nt+1, ' of ', len(THETA)
			FD.n += 1

		MatrizT[:,nt]	= np.copy(FD.receiver_signals[:,1])
		MatrizR[:,nt]	= np.copy(FD.receiver_signals[:,0])
		nt += 1
	
	stop    = time.time()
	print stop-start, " seconds"
	
	data={}
	data['Transmission'] = MatrizT
	data['Reflection']  = MatrizR
	
	# Save in Matlab Files (.mat)
	savemat(file_data,data)
	
	figure(1); imshow(MatrizT,aspect='auto'); 
	show()
	
	
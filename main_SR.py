"""
Execution for a single angle rotation of the sample

Miguel Molero (M. Molero)  & Ursula Iturraran-Viveros, 2012

miguel.molero@gmail,com
ursula.iturraran@gmail.com

"""

from   ANI2D_Classes import *
from   ANI2D         import *
import sys
import time
from   scipy.io      import savemat

file_data = 'test_0'

Scale  = 1.0

Width  = Scale*(75)
Height = Scale*(100)
Tap    = 6

Width_Sample   = Scale*(35)   # Dimension in mm
Height_Sample  = Scale*(50)   # Dimension in mm

Frequency      = 500e3   # Operating Frequency
#Simtime is not applied, time step is set using  Ntime
SimTime        = 100e-6; # Simulaton Time
#Simulation step time
Ntime          = 5000

device         = 'GPU'    # 'CPU' or 'CPU'

THETA = [0]    #for a single rotation angle


if __name__ == '__main__':

	from  matplotlib.pylab   import *

	#Create Scenario
	Scenario    = NewImage(Width,Height,Pixel_mm = 10, GrayScale=0)
	Scenario.I  = addRectangle(Scenario, int(Width*0.5), int(Height*0.5), Width_Sample, Height_Sample, 60, 0)
	
	#Material Setup
	M1         = Material( name="M1",rho=1000,c11=2.19e9,c12=2.19e9,c22=2.19e9,c44=1e-30,gray=0)
	M2         = Material( name="M2",rho=2648,c11=15.9e9,c12=7.0e9,c22=15.9e9,c44=3.4e9,gray= 60)
	materials = MaterialSet()
	materials.append(M1)
	materials.append(M2)
		
	#Absorbing conditions
	Scenario.createABC(Tap)
	
	#Create Inspection Type, Dual = Reflection & Transmission Inspection
	source        = Source(LaunchType = 'Dual')
	emitter       = Transducer(Size=50,  Offset=0)
	receiver      = Transducer(Size=2.0, Offset=0)
	
	transducer    = [emitter, receiver]
	
	waveform      = Signal(Amplitude=1, Frequency=Frequency)
	
	
	sim_model    = SimulationModel(TimeScale=1, MaxFreq=2*Frequency, PointCycle=10, SimTime=SimTime)
	sim_model.jobParameters(materials)
	sim_model.createNumericalModel(Scenario)
	sim_model.initReceivers()
	sim_model.setDevice(device)
	
	print sim_model.MRI, sim_model.NRI
	
	#########
	###  Ntime x sim_model.Ntime
	#Transmission and Reflection Matrices 
	MatrizT = np.zeros((Ntime,np.size(THETA,0)),dtype=np.float32)
	MatrizR = np.zeros((Ntime,np.size(THETA,0)),dtype=np.float32)
	#########
	
	
	nt = 0
	
	#Create Simulation Object
	FD = ANI2D(Scenario, materials, source, transducer, waveform, sim_model)
	
	##########
	##########
	#Array used to control iteration number
	t = sim_model.dt*np.arange(0,Ntime)
	#Init waveform for emitter transducer
	FD.Init_Source(waveform, t)
	FD.AgainReceiverSetup(Scenario,source,transducer, sim_model,Ntime)
	#########

	start = time.time()
	while (FD.n < Ntime):
		FD.RunSim()
		if FD.n % 500==0:
			print FD.n, " of iteration: ", Ntime
		FD.n += 1

	stop    = time.time()
	print stop-start, " seconds"

	MatrizT[:,nt]	= np.copy(FD.receiver_signals[:,1])
	MatrizR[:,nt]	= np.copy(FD.receiver_signals[:,0])
	nt += 1
	

	data={}
	data['Transmission'] = MatrizT
	data['Reflection']  = MatrizR
	
	# Save in Matlab Files (.mat)
	savemat(file_data,data)
	

	figure(1); plot(MatrizT); 
	show()
	
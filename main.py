"""Ursula Iturraran-Viveros  & Miguel Molero, 2013

   ursula.iturraran@gmail.com, miguel.molero@gmail.com
"""

version = "0.11"

import	 sys
import	 time
import	 numpy	as	 np

import	 matplotlib.pyplot	 as plt
import   glumpy


Plotting    = True

EnableVideo = True #only works when plotting
FILE        = "test"  # FILE name, save in .mat format



if __name__ == '__main__':

	from	 FD25D_CL	 import *

	# define scenario dimension limits
	scenario  = [3.0, 0, 8] # rmax, zmin, zmax

	# define source parameters and source position: n=0 => Monopolar, n=1 =>Dipolar
	source    = [0, 8000.0, 0.005, 4] 	# n, fmax, r, z

	#define media parameters

#       This case corresponds to the Fast formation Fig. 3. Uncomment to use it an comment the 
#       following set of parameters for Fig. 5
#	rx  =	[0.1,    3.0]   #radius borehole
#	vpx = [1500.0, 4000.0]  #longitudinal velocity
#	vsx = [1e-30,  2300.0]  #shear velocity
#	rho = [1000.0, 2500.0]  #density

#       This case corresponds to Fig. 5 where there is an invaded zone
	rx  =	[0.1,  0.5,  3.0]   #radius borehole and for the rings
	vpx = [1500.0, 3000, 4000.0]  #longitudinal velocity
	vsx = [1e-30,  1800, 2300.0]  #shear velocity
	rho = [1000.0, 2300, 2500.0]  #density
        
	# define Time Iteration
	TimeIter	= 6000
	# absorbing layer size (absorbing boundaries)
	Tap			= 20

	# Define Position of Receivers
	pz = np.linspace(6.75,7.8,8,endpoint=True)
	pr = [0.01]*np.size(pz,0)

	########

	FD = FD25D(scenario, source, rx, vpx, vsx, rho)
	FD.InitCL("GPU")  # change for "CPU" when running on a CPU
	FD.Setup(Tap)
	####if not interested the Layer intercepting just FD.MediumSetup()
	FD.MediumSetup(LayerInterProps=[2200, 1800, 967],LayerInterDim=[3.9, 4.2],LayerIntercepting=True)
	FD.absSetup()
	FD.receivers_setup(pr, pz, TimeIter)
	FD.Init_Fields_CL()


	start = time.time()

  	if Plotting:
  
		Z		 =	FD.SV.transpose()
		fig		 =	glumpy.figure((int(FD.NRI*0.5),int(FD.NZI*0.5)) )
		I		 =	glumpy.Image(Z, interpolation='bilinear', colormap= glumpy.colormap.IceAndFire,
								vmin=-50, vmax=0)
																	
		@fig.event
		def on_key_press(key, modifiers):
			if key == glumpy.window.key.ESCAPE:
				sys.exit();
			else:
				pass

		@fig.event
		def on_draw():
			fig.clear()
			I.draw(x=0, y=0, z=0, width=fig.window.width, height=fig.window.height)


		@fig.event
		def on_idle(*args):

			if (FD.t < TimeIter):
				FD.RunCL()
			
				if FD.t % 200==0:
					print str(FD.t)  + " Total " + str(TimeIter)

				if (FD.t %100==0):
					FD.GL()
					Z[...]	=  FD.SV.transpose()
					I.update()
					fig.redraw()
				
					if EnableVideo:
						fileName = FILE + str(FD.t/100) + ".jpg"
						FD.save_video(fig,fileName)

			else:
					#Retrieve the receiver signals from the computing device to the host
					FD.saveOutput() #FD.receivers_signals
					# Save simulation data: receivers, dt, dr, dz
					FD.save_data(FILE)
					sys.exit();


		glumpy.show()






	else:
		# main loop
		while (FD.t < TimeIter):
			FD.RunCL()
			if FD.t % 500 == 0:
				print FD.t, " of total iterations: ", TimeIter
		print time.time()-start


		#Retrieve the receiver signals from the computing device to the host
		FD.saveOutput() #FD.receivers_signals

		# Save simulation data: receivers, dt, dr, dz
		FD.save_data(FILE)

		# plot receivers
		plt.figure();
		ind = 0
		sig = []

		vector_t = np.arange(0,TimeIter)*FD.dt

		N	  = np.size(FD.receivers_signals,1)
		delta = 2.0/float(N)

		for i in range(0,N):
			sig = FD.receivers_signals[:,i]/np.max(FD.receivers_signals)
			plt.plot(vector_t,sig+ind)
			ind += delta
			plt.hold(True)
		plt.show()



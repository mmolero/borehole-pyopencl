"""Ursula Iturraran  & Miguel Molero, 2012

   ursula.iturraran@gmail.com, miguel.molero@gmail.com
"""

version = "0.11"

import	 numpy	as	         np
from	 math	import       pi, exp

from     scipy.weave         import inline
from     scipy.io            import savemat

import   pyopencl            as cl
import   glumpy

from     PIL       import Image
import	 OpenGL.GL as     GL

 

def	  ricker(t,ts,fsavg):
	# ricker pulse
	a  = fsavg*pi*(t-ts)
	a2 = a*a
	return (1.0-2.0*a2)*exp(-a2)


class FD25D:

	def __init__(self, scenario, source,  rx, vpx, vsx, rho):

		self.rmax,self.zmin,self.zmax                = scenario
		self.n, self.fsavg,self.rsource,self.zsource = source

		self.n = np.float32(self.n)

		self.nrings = len(rx)
		self.rx	  	= np.array(rx)
		self.vpx  	= np.array(vpx)
		self.vsx  	= np.array(vsx)
		self.rho 	= np.array(rho)

		self.t = 0

	def InitCL(self, DEVICE="CPU"):

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


	def Setup(self,Tap):

		self.Tap = Tap

		#Determine maximum temporal frequency to avoid dispersion
		#Determine mininum spatial sampling interval
		#10 empirical sample per wavelength	 fmax=favg*4

		if (self.vsx[1])<1:
			vmin = np.min(self.vpx)
		else:
			vmin = np.min([ np.min(self.vpx), self.vsx[1] ] )


		vmax = np.max([ self.vpx ])
		
		
		#Determine spatial grid
		h  = vmin / ( 10.0 * self.fsavg * 4.)
		print "dr, dz = ", h
		self.dr = np.float32(h)
		self.dz = np.float32(h)

		h = min([np.abs(self.dr), np.abs(self.dz) ])

		#Determine time sampling interval to ensure stability, this is cartesian stablity
		self.dt	  = 0.25 * h / ( vmax )
		print "dt = ", self.dt

		self.NR =  int (self.rmax/self.dr )
		self.NZ =  int((self.zmax - self.zmin)/self.dz)

		self.i0	  = int(self.rsource/self.dr)
		self.j0	  = int((self.zsource-self.zmin)/self.dr)


		if	self.i0 < 0:
			self.i0 = 0

		#Dipole should not be at zero
		if( self.i0 < 1 and	 self.n > 0):
			self.i0 = 1
		if( self.j0 < 0):
			self.j0 = 0
		if( self.j0 > self.NZ-6):
			self.j0 = self.NZ-6

		#Parameters for the Ricker wavelet
		self.ts = 1.0 / self.fsavg

	def MediumSetup(self, LayerInterProps,LayerInterDim,LayerIntercepting=False):
		#Medium parameters

		self.rx[self.nrings-1] = (self.NR) * self.dr
		#From outside in and overwrite properties

		rz = self.dz*np.arange(0,self.NZ)
		rr = self.dr*np.arange(0,self.NR)

		numrings = range(0,self.nrings)

		self.den_tmp	= np.zeros((self.NR,self.NZ), dtype=np.float32)
		self.mu_tmp		= np.zeros((self.NR,self.NZ), dtype=np.float32)
		self.lam_tmp	= np.zeros((self.NR,self.NZ), dtype=np.float32)

		for k in reversed(numrings):

			indr = np.nonzero(rr<=self.rx[k])
			self.den_tmp[indr,:]  = self.rho[k]
			self.mu_tmp [indr,:]  = self.rho[k]*(self.vsx[k]**2)
			self.lam_tmp[indr,:]  = self.rho[k]*(self.vpx[k]**2) - 2.*self.rho[k]*(self.vsx[k]**2)

        ########
		if LayerIntercepting:
			indr = np.ceil(self.rx[0]/self.dr);
			print indr
			z1 = LayerInterDim[0]
			z2 = LayerInterDim[1]
			for i in range(0, self.NZ):
				if ( (i*self.dz > z1) and (i*self.dz < z2) ):
					self.den_tmp[indr:,i] = LayerInterProps[0]
					self.mu_tmp [indr:,i] = LayerInterProps[0]*(LayerInterProps[2]**2)
					self.lam_tmp[indr:,i] = LayerInterProps[0]*(LayerInterProps[1]**2) - 2.*LayerInterProps[0]*(LayerInterProps[2]**2)
      	########

			

		self.NRI = self.NR + self.Tap
		self.NZI = self.NZ + 2*self.Tap
		print "NRI = ", self.NRI, " NZI = ", self.NZI
		print "i0 = ", self.i0, "  j0 = ", self.j0

		self.den	= np.ones((self.NRI,self.NZI), dtype=np.float32)*np.float32(self.rho[-1])
		self.mu		= np.ones((self.NRI,self.NZI), dtype=np.float32)*np.float32(self.rho[-1]*(self.vsx[-1]**2))
		self.lam	= np.ones((self.NRI,self.NZI), dtype=np.float32)*np.float32(self.rho[-1]*(self.vpx[-1]**2) - 2.*self.rho[-1]*(self.vsx[-1]**2))

		self.den[0 : self.NRI-self.Tap, self.Tap : self.NZI-self.Tap]  = np.copy(self.den_tmp)
		self.mu [0 : self.NRI-self.Tap, self.Tap : self.NZI-self.Tap]  = np.copy(self.mu_tmp)
		self.lam[0 : self.NRI-self.Tap, self.Tap : self.NZI-self.Tap]  = np.copy(self.lam_tmp)
		

		#Arrays of stresses
		self.srr = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.szz = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.stt = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.srz = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.stz = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.srt = np.zeros((self.NRI,self.NZI), dtype=np.float32)

		#Arrays of velocities
		self.vr = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.vt = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.vz = np.zeros((self.NRI,self.NZI), dtype=np.float32)

		#Elastic constants
		self.lam2mu = np.zeros((self.NRI,self.NZI), dtype=np.float32)

		#Densities
		self.denr	= np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.denz	= np.zeros((self.NRI,self.NZI), dtype=np.float32)

		#Shear modulus harmonic average
		self.mur	 = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.muz	 = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.murz	 = np.zeros((self.NRI,self.NZI), dtype=np.float32)

		self.ABSC	 = np.zeros((self.NRI,self.NZI), dtype=np.float32)
		self.SV		 = np.zeros((self.NRI,self.NZI), dtype=np.float32)

		#Harmonic average for mu=mur,muz, average for density
		self.mur  = np.copy(self.mu)
		self.muz  = np.copy(self.mu)
		self.murz = np.copy(self.mu)
		self.denr = np.copy(self.den)
		self.denz = np.copy(self.den)

		self.denr[1:,:]= ( self.den[1:,:] + self.den[0:-1,:]  )/2.
		self.denz[:,1:]= ( self.den[:,1:] + self.den[:,0:-1]  )/2.

		mu   = np.copy(self.mu);
		mur  = np.copy(self.mur);
		muz  = np.copy(self.muz);
		murz = np.copy(self.murz);

		NZI = self.NZI
		NRI = self.NRI

		extra_code = """
				       #define ind(i,j,NZI)     ( ((i)*(NZI))+ (j) )

				"""

		code  = """

			for(int j=1; j<NZI; ++j){
				for(int i=1; i<NRI; ++i){

					if ( ( mu[ind(i,j,NZI)] == 0.0 ) || (mu[ind(i-1,j,NZI)] == 0.0) ){
						mur[ind(i,j,NZI)] = 0.0;
					}

					else{
						mur[ind(i,j,NZI)] = 2.0*mu[ind(i,j,NZI)]*mu[ind(i-1,j,NZI)]/(mu[ind(i,j,NZI)]+mu[ind(i-1,j,NZI)] );
					}

					if ( ( mu[ind(i,j,NZI)] == 0.0 )  || (mu[ind(i,j-1,NZI)]==0.0) ){
						muz[ind(i,j,NZI)] = 0.0;
					}
					else{
						muz[ind(i,j,NZI)] = 2.0*mu[ind(i,j,NZI)]*mu[ind(i,j-1,NZI)]/( mu[ind(i,j,NZI)] + mu[ind(i,j-1,NZI)] );
				    }

				   if ( (murz[ind(i,j,NZI)] ==0.0) ||   (mu[ind(i-1,j,NZI)]== 0.0) || (mu[ind(i,j-1,NZI)] == 0.0) || (mu[ind(i-1,j-1,NZI)]==0.0))
				   {
				     murz[ind(i,j,NZI)] = 0.0 ;
				   }

			       else{

						double tmp1 =  mu[ind(i,j,NZI)]   * mu[ind(i-1,j,NZI)] ;
						double tmp2 =  mu[ind(i,j-1,NZI)] * mu[ind(i-1,j-1,NZI)];
						double tmp3 =  (mu[ind(i,j,NZI)]  * mu[ind(i-1,j,NZI)] * mu[ind(i,j-1,NZI)]);
						double tmp4 =  (mu[ind(i,j,NZI)]  * mu[ind(i-1,j,NZI)] * mu[ind(i-1,j-1,NZI)]);
						double tmp5 =  (mu[ind(i,j,NZI)]  * mu[ind(i,j-1,NZI)] * mu[ind(i-1,j-1,NZI)]);
						double tmp6 =  (mu[ind(i-1,j,NZI)]* mu[ind(i,j-1,NZI)] * mu[ind(i-1,j-1,NZI)]);

						double tmp7    = 4.0*tmp1*tmp2;
						double tmp8    = tmp3 + tmp4 + tmp5 + tmp6;
						murz[ind(i,j,NZI)] = tmp7/tmp8;
				    }
				  }
			     }

				"""

		inline(code,['mu', 'mur', 'muz', 'murz','NRI','NZI'], support_code = extra_code, compiler='gcc')

		self.mur  = np.copy(mur);
		self.muz  = np.copy(muz);
		self.murz = np.copy(murz);

		#Take the inverse of the average densities, multiply elastic constant by dt and make lam2mu
		self.den   = np.float32(self.dt/self.den)
		self.denr  = np.float32(self.dt/self.denr)
		self.denz  = np.float32(self.dt/self.denz)

		self.lam2mu[:,:]  = np.float32(self.dt*( self.lam[:,:] + 2.0*self.mu[:,:] ))
		self.lam[:,:]	  = np.float32(self.dt*self.lam[:,:])
		self.mur[:,:]	  = np.float32(self.dt*self.mur[:,:])
		self.muz[:,:]	  = np.float32(self.dt*self.muz[:,:])
		self.murz[:,:]	  = np.float32(self.dt*self.murz[:,:])

		#Arrays for radius
		self.rii	  = np.zeros((1,self.NRI), dtype=np.float32)
		self.rhalf	  = np.zeros((1,self.NRI), dtype=np.float32)

		#Take r an array and avoid divides
		for i in range(1, self.NRI-self.Tap):
			self.rii[0,i]   = 1.0/( (i)*self.dr )

		for i in range(0, self.NRI-self.Tap):
			self.rhalf[0,i] = 1.0/( (i-0.5)*self.dr )

		#Take the inverse of dz to avoid divides
		self.drii = np.float32(1.0/self.dr)
		self.dzjj = np.float32(1.0/self.dz)


	def absSetup(self):
		#ABSORBING	BOUNDARY CONDITION	see Cerjan (1985)
		#Note that the left hand side is not an ABC because it is where
		#he borehole axis is located
		APARA = 0.015
		for i in range(0,self.NRI):
			for j in range(0,self.NZI):

				if	 i	< self.Tap:
					self.ABSC[i,j] = 1.0
				elif j < self.Tap:
					self.ABSC[i,j] = np.exp( -( ( APARA * (self.Tap-j))**2 ) )
				elif  i >  self.NRI-self.Tap+1:
					self.ABSC[i,j] = np.exp( -( ( APARA * (i-self.NRI+self.Tap-1) )**2 ) )
				elif  j > self.NZI-self.Tap+1:
					self.ABSC[i,j] = np.exp( -( ( APARA * (j-self.NZI+self.Tap-1) )**2 ) )
				else:
					self.ABSC[i,j] = 1.0



	def save_data(self, File):

		data = {}
		data['receivers'] = np.copy(self.receivers_signals)
		data['dt']        = np.copy(self.dt)
		data['dr']        = np.copy(self.dr)
		data['dz']        = np.copy(self.dz)

		savemat(File,data)
		
	
	def save_video(self, fig, File):
		# Capture image from the OpenGL buffer
		buffer = (GL.GLubyte * (3*fig.window.width*fig.window.height) )(0)
		GL.glReadPixels(0, 0, fig.window.width, fig.window.height, GL.GL_RGB, GL.GL_UNSIGNED_BYTE, buffer)

		# Use PIL to convert raw RGB buffer and flip the right way up
		image = Image.fromstring(mode="RGB", size=(fig.window.width, fig.window.height), data=buffer)	   
		image = image.transpose(Image.FLIP_TOP_BOTTOM)
		image.save(File)
		
		

	def Init_Fields_CL(self):

		self.srr_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.srr)
		self.stt_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.stt)
		self.szz_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.szz)
		self.srt_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.srt)
		self.srz_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.srz)
		self.stz_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.stz)

		self.vr_buf		    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.vr)
		self.vt_buf		    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.vt)
		self.vz_buf		    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.vz)

		self.den_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.den)
		self.denr_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.denr)
		self.denz_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.denz)

		self.rii_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.rii)
		self.rhalf_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.rhalf)
		self.ABSC_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.ABSC)

		self.lam2mu_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.lam2mu)
		self.lam_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.lam)
		self.murz_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.murz)
		self.mur_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.mur)
		self.muz_buf	    = cl.Buffer(self.ctx, self.mf.READ_WRITE	 | self.mf.COPY_HOST_PTR, hostbuf=self.muz)

		self.receivers_buf  = cl.Buffer(self.ctx, self.mf.READ_WRITE | self.mf.COPY_HOST_PTR, hostbuf=self.receivers_signals)

		self.receiver_r_buf = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=self.receiver_r)
		self.receiver_z_buf = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=self.receiver_z)

		self.program	    = cl.Program( self.ctx, self.Kernel_CL() ).build()

		print "n = ", self.n


	def receivers_setup(self, pr, pz, Tiempo):

		self.receiver_r         = np.int32(np.around(np.array(pr)/self.dr))
		self.receiver_z         = np.int32(np.around(np.array(pz)/self.dz))

		self.N_z                = np.size(self.receiver_z,0)
		self.receivers_signals  = np.zeros((Tiempo,self.N_z), dtype=np.float32)


	def get_receivers(self,t):

		cl.enqueue_copy(self.queue, self.srr, self.srr_buf)

		for cxy in range(0,self.N_z):
			self.receivers_signals[t,cxy] = self.srr[self.receiver_r[cxy],self.receiver_z[cxy]]


	def RunCL(self):
          #Source
          ti     = self.t*self.dt
          fuente = np.float32(ricker(ti,self.ts,self.fsavg))

          #computing velocities
          self.program.velocity_FD25D(  self.queue, (self.NZI,self.NRI,), None,
											  self.srr_buf,
											  self.stt_buf,
											  self.szz_buf,
											  self.stz_buf,
											  self.srz_buf,
											  self.srt_buf,
											  self.vr_buf,
											  self.vt_buf,
											  self.vz_buf,
											  self.den_buf,
											  self.denr_buf,
											  self.denz_buf,
											  self.rii_buf,
											  self.rhalf_buf,
											  self.ABSC_buf,
											  fuente).wait()
          #computing stresses
          self.program.stress_FD25D(self.queue, (self.NZI,self.NRI,), None,
												 self.srr_buf,
												 self.stt_buf,
												 self.szz_buf,
												 self.stz_buf,
												 self.srz_buf,
												 self.srt_buf,
												 self.vr_buf,
												 self.vt_buf,
												 self.vz_buf,
												 self.lam2mu_buf,
												 self.lam_buf,
												 self.mur_buf,
												 self.muz_buf,
												 self.murz_buf,
												 self.rii_buf,
												 self.rhalf_buf,
												 self.ABSC_buf,
                                                 self.receivers_buf,
								                 self.receiver_r_buf,
								                 self.receiver_z_buf,
								                 np.int32(self.t) ).wait()
          self.t +=1

		 
		 
	def GL(self):
	
		cl.enqueue_copy(self.queue, self.srr, self.srr_buf)
		self.SV =  self.srr
		self.SV = 20.*np.log10((np.abs(self.SV)/np.max(np.abs(self.SV+1e-40))) + 1e-40)



	def saveOutput(self):
		cl.enqueue_copy(self.queue, self.receivers_signals, self.receivers_buf)

	def Kernel_CL(self):

		extra_code	 =	"""

					 #define ind(i,j)  ( ((NZI)*(i)) + (j) )
					 #define NRI	   %s
					 #define NZI	   %s
					 #define n		   %g
					 #define drii	   %g
					 #define dzjj	   %g
					 #define i0		   %s
					 #define j0		   %s
					 #define Nz        %s
					 #define ir(i,j)   (  ( (Nz)*(i) ) + (j)     )
		"""%( str(self.NRI),str(self.NZI), np.float32(self.n), np.float32(self.drii), np.float32(self.dzjj),
		      str(self.i0), str(self.j0), str(self.N_z) )



		kernel_source = extra_code + r"""


				__kernel void  get_receivers( __global float *srr,
											  __global float *receivers,
				     						        __global const int *receiver_r,
											  __global const int *receiver_z,
											  int t){

					uint m,r,z;
					m  =  get_global_id(0);
					r  = receiver_r[m];
					z  = receiver_z[m];
					receivers[ir(t,m)] = srr[ind(r,z)];
				}


				 __kernel void stress_FD25D(__global float *srr,
											__global float *stt,
											__global float *szz,
											__global float *stz,
											__global float *srz,
											__global float *srt,
											__global float *vr,
											__global float *vt,
											__global float *vz,
											__global  float *lam2mu,
											__global  float *lam,
											__global  float *mur,
											__global  float *muz,
											__global  float *murz,
											__global  float *rii,
											__global  float *rhalf,
											__global  float *absc,
                                            __global  float  *receivers,
								            __global  const int    *receiver_r,
								            __global  const int    *receiver_z,
								             int  t){



					 int j = get_global_id(0);
				     int i = get_global_id(1);

                            uint m,r,z;


					  if (n == 0.0){

						    if (j>0 && j <NZI-1){

						      srr[ind(0,j)] +=	(   (lam2mu[ind(0,j)] * drii)	        *   vr[ind(1,j)]  +
											        (lam   [ind(0,j)] * dzjj)	        * ( vz[ind(0,j+1)]  - vz[ind(0,j)] )+
										            (lam   [ind(0,j)] * rhalf[0]*0.5f)  *   vr[ind(1,j)] 	 +
											        (lam   [ind(0,j)] * rhalf[0])       *  ( n * vt[ind(0,j)] ) );


						      stt[ind(0,j)] +=	(   (lam   [ind(0,j)] * drii)	        *   vr[ind(1,j)]	 +
											        (lam   [ind(0,j)] * dzjj)	        * ( vz[ind(0,j+1)]  - vz [ind(0,j)] ) +
										            (lam2mu[ind(0,j)] * 0.5f*rhalf[0])  *   vr[ind(1,j)]	 +
											        (lam2mu[ind(0,j)] * rhalf[0] )      * ( n * vt[ind(0,j)] ) );

						      szz[ind(0,j)] +=	(   (lam   [ind(0,j)] * drii)	        *   vr[ind(1,j)]  +
											        (lam2mu[ind(0,j)] * dzjj)	        * ( vz[ind(0,j+1)] -	vz[ind(0,j)] ) +
										            (lam   [ind(0,j)] * 0.5f *rhalf[0]) *   vr[ind(1,j)]	 +
											        (lam   [ind(0,j)] * rhalf[0])       * ( n * vt[ind(0,j)] ) );
                           } // end if
					  }// end if

					  else{

							if ( j > 0 ){

							     srr[ind(0,j)] = 0.0f;
							     stt[ind(0,j)] = 0.0f;
							     szz[ind(0,j)] = 0.0f;

							     stz[ind(0,j)] += ( (muz [ind(0,j)] * dzjj)  * ( vt[ind(0,j)]	 - vt[ind(0,j-1)] )	 -
							                        (muz [ind(0,j)] * rii[i] * n  * vz[ind(0,j)] ) );

							     srz[ind(0,j)] += ( (murz[ind(0,j)]* drii)   * ( vz[ind(0,j)] )  +
							                        (murz[ind(0,j)]* dzjj)   * ( vr[ind(0,j)] - vr[ind(0,j-1)] ) );

							     srt[ind(0,j)] = 0.0f;

							} // end if
					   } // end else


					barrier(CLK_GLOBAL_MEM_FENCE);


					 if(i>0 && i< NRI-1 && j> 0 && j<NZI-1) {

								srr[ind(i,j)]	  +=	( (lam2mu [ind(i,j)]  *  drii )          * (  vr[ind(i+1,j)] - vr[ind(i,j)]  )  +
														  (lam	  [ind(i,j)]  *  dzjj )          * (  vz[ind(i,j+1)] - vz[ind(i,j)]  )  +
														  (lam	  [ind(i,j)]  *  rhalf[i]*0.5f ) * (  vr[ind(i+1,j)] - vr[ind(i,j)]  )  +
														  (lam	  [ind(i,j)]  *  rhalf[i] )	     * (  n * vt[ind(i,j)] )  );



								stt[ind(i,j)]	  +=	( (lam   [ind(i,j)]   *  drii )           * (  vr[ind(i+1,j)] - vr[ind(i,j)]  ) +
														  (lam   [ind(i,j)]   *  dzjj )           * (  vz[ind(i,j+1)] - vz[ind(i,j)]  ) +
														  (lam2mu[ind(i,j)]   * rhalf[i] * 0.5f)  * (  vr[ind(i+1,j)] - vr[ind(i,j)]  ) +
														  (lam2mu[ind(i,j)]   * rhalf[i] )        * (  n * vt[ind(i,j)] ) );

								szz[ind(i,j)]	  +=	( (lam	 [ind(i,j)]   * drii )	          * (  vr[ind(i+1,j)] - vr[ind(i,j)]  ) +
														  (lam2mu[ind(i,j)]   * dzjj )	          * (  vz[ind(i,j+1)] - vz[ind(i,j)]  ) +
														  (lam	 [ind(i,j)]   * rhalf[i] * 0.5f ) * (  vr[ind(i+1,j)] - vr[ind(i,j)]  ) +
														  (lam	 [ind(i,j)]   * rhalf[i] )		  * (  n * vt[ind(i,j)] ) );

					 }


					 if(i>0 && i< NRI-1 && j> 0 && j<NZI-1) {

								stz[ind(i,j)]	 += (	 (muz [ind(i,j)] *  dzjj)          * (  vt[ind(i,j)]	- vt[ind(i,j-1)] ) -
														 (muz [ind(i,j)] * rii [i] * n)	   *	vz[ind(i,j)]  );

								srz[ind(i,j)]	 += (	 (murz[ind(i,j)] *  drii)          * (  vz[ind(i,j)]	- vz[ind(i-1,j)] ) +
														 (murz[ind(i,j)] *  dzjj)          * (  vr[ind(i,j)]    - vr[ind(i,j-1)] ) );

								srt[ind(i,j)]	 += (	 (mur [ind(i,j)] * drii)           * (  vt[ind(i,j)]	- vt[ind(i-1,j)] ) -
													     (mur [ind(i,j)] * rii [i] * 0.5f) * (  vt[ind(i,j)]	- vt[ind(i-1,j)] )-
														 (mur [ind(i,j)] * rii [i] *	n) * (  vr[ind(i,j)])  );

					 }


				  barrier(CLK_GLOBAL_MEM_FENCE);

					   srr[ind(i,j)] *= absc[ind(i,j)];
					   stt[ind(i,j)] *= absc[ind(i,j)];
					   szz[ind(i,j)] *= absc[ind(i,j)];
					   srz[ind(i,j)] *= absc[ind(i,j)];
					   srt[ind(i,j)] *= absc[ind(i,j)];
					   stz[ind(i,j)] *= absc[ind(i,j)];

                         barrier(CLK_GLOBAL_MEM_FENCE);

                        for ( m=0; m<Nz; m++){
                             r  = receiver_r[m];
					z  = receiver_z[m];
					receivers[ir(t,m)] = srr[ind(r,z)];
                         }


			   }


	     __kernel void velocity_FD25D(__global float *srr,
								      __global float *stt,
								      __global float *szz,
								      __global float *stz,
								      __global float *srz,
								      __global float *srt,
								      __global float *vr,
								      __global float *vt,
								      __global float *vz,
								      __global  float *den,
							            __global  float *denr,
								      __global  float *denz,
								      __global  float *rii,
								      __global  float *rhalf,
								      __global  float *absc,
								        float source){


					int j = get_global_id(0);
				    int i = get_global_id(1);


					 if ( n == 0  ){

							if ( j > 0 ){
									vz[ind(0,j)] += ( (denz[ind(0,j)] * 0.5f* rhalf[0]) * ( srz[ind(0,j)]	 + 2.f * n * stz[ind(0,j)] ) +
													   denz[ind(0,j)] * drii            *   srz[ind(0,j)]    +
						     						   denz[ind(0,j)] * dzjj            * ( szz[ind(0,j)] - szz[ind(0,j-1)] ) );
						 	}
					  }


					 else {
							 if (j > 0 && j<NZI-1){

								vr[ind(0,j)] +=	 ( (denr[ind(0,j)] *rii[i]) * ( srr[ind(1,j)] - stt[ind(1,j)] ) +
												   (denr[ind(0,j)] *drii  ) * ( srr[ind(1,j)] ) +
												   (denr[ind(0,j)] *dzjj  ) * (	srz[ind(0,j+1)]- srz[ind(0,j)] ) );

								vt[ind(0,j)] += (  (den[ind(0,j)]  *rhalf[0] ) *  srt[ind(1,j)]  +
								                   (denr[ind(0,j)] *drii     ) *  srt[ind(1,j)]  +
												   (denr[ind(0,j)] *dzjj     ) *( stz[ind(0,j)] - stz[ind(0,j-1)] ) );

								vz[ind(0,j)] = 0.0f;
							 }
					  }


				    barrier(CLK_GLOBAL_MEM_FENCE);




					if( i > 1 && i< NRI-1 && j > 1 && j<NZI-1 ){


					   	  vr[ind(i,j)] += ( (denr[ind(i,j)] * 0.5f*rii[i]) *( srr[ind(i,j)]   + srr[ind(i-1,j)] - stt[ind(i,j)] - stt[ind(i-1,j)]  +  2.f* n * srt[ind(i,j)] ) +
					                        (denr[ind(i,j)] * drii)        *( srr[ind(i,j)]   - srr[ind(i-1,j)]) +
					                        (denr[ind(i,j)] * dzjj)        *( srz[ind(i,j+1)] - srz[ind(i,j)]) );


					      vt[ind(i,j)] +=  ( (den[ind(i,j)] * rhalf[i])    * ( srt[ind(i+1,j)] + srt[ind(i,j)] - n * stt[ind(i,j)] ) +
										   (  den[ind(i,j)] * drii)        * ( srt[ind(i+1,j)] - srt[ind(i,j)] )  +
										   (  den[ind(i,j)] * dzjj)        * ( stz[ind(i,j+1)] - stz[ind(i,j)] ) );

  						  vz[ind(i,j)] +=  ( (denz[ind(i,j)] *  0.5f*rhalf[i])*( srz[ind(i+1,j)] + srz[ind(i,j)] + 2.0f* n * stz[ind(i,j)] ) +
											 (denz[ind(i,j)] *	drii)*( srz[ind(i+1,j)] - srz[ind(i,j)] ) +
											 (denz[ind(i,j)] *  dzjj)*( szz[ind(i,j)]   - szz[ind(i,j-1)] )  );

					}


					 barrier(CLK_GLOBAL_MEM_FENCE);


					    vr[ind(i,j)] *= absc[ind(i,j)];
						vt[ind(i,j)] *= absc[ind(i,j)];
						vz[ind(i,j)] *= absc[ind(i,j)];

					if (n==0){
					    srr[ind(i0,j0)] -= source;
					    stt[ind(i0,j0)] -= source;
					    szz[ind(i0,j0)] -= source;
					}
					else{
					     vr[ind(i0,j0)] -= source;

					}


		    }

		"""

		return kernel_source








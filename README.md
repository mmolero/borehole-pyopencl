###Borehole-pyopencl

if you use borehole-pyopencl code in your research, we would appreciate the citation of the following article:

**"Simulation of sonic waves along a borehole in heterogeneous formation: Accelerating 2.5-D finite differences using PyOpenCL"**, Ursula Iturrarán-Viveros, Miguel Molero, Computer & Geosciences 56, 2013, 161-169

[link](http://www.sciencedirect.com/science/article/pii/S0098300413000782) 

____

Abstract:

This article presents an implementation of a 2.5-D finite- difference (FD) code to model acoustic full waveform monopole logging in cylindrical coordinates accelerated by using the new parallel computing devices (PCDs). For that purpose we use the industry open standard Open Computing Language (OpenCL) and an open-source toolkit called PyOpenCL. The advantage of OpenCL over similar languages is that it allows one to program a CPU (Central Processing Unit) a GPU (Graphics Processing Unit), or multiple GPUs and their interaction among them and with the CPU, or host device. We describe the code and give a performance test in terms of speed using six different computing devices under different operating systems. A maximum speedup factor over 34.2, using the GPU is attained when compared with the execution of the same program in parallel using a CPU quad-core. Furthermore, the results obtained with the finite differences are validated using the discrete wavenumber method (DWN) achieving a good agreement. To provide the Geoscience and the Petroleum Science communities with an open tool for numerical simulation of full waveform sonic logs that runs on the PCDs, the full implementation of the 2.5-D finite difference with PyOpenCL is included. 

____

Some examples:

  - Monopolar source and a fast formation around the fluid-filled borehole
  
  <a href="http://www.youtube.com/watch?feature=player_embedded&v=ieVFxDvveDQ
  " target="_blank"><img src="http://img.youtube.com/vi/ieVFxDvveDQ/0.jpg" 
  alt="Monopolar source and a fast formation around the fluid-filled borehole"    
  width="240" height="240" border="10" /></a>  

  - Monopolar source and a slow formation around the fluid-filled borehole
  
  <a href="http://www.youtube.com/watch?feature=player_embedded&v=8OMLNzg79sI
  " target="_blank"><img src="http://img.youtube.com/vi/8OMLNzg79sI/0.jpg" 
  alt="Monopolar source and a slow formation around the fluid-filled borehole"    
  width="240" height="240" border="10" /></a>

  - Monopolar source and a fast formation around the fluid-filled borehole and the borehole is intercepted by a slow       layer 
  
  <a href="http://www.youtube.com/watch?feature=player_embedded&v=HxPb6lMvvpY
  " target="_blank"><img src="http://img.youtube.com/vi/HxPb6lMvvpY/0.jpg" 
  alt="Monopolar source and a fast formation around the fluid-filled borehole and the borehole is intercepted by a slow    layer"    
  width="240" height="240" border="10" /></a>
  

___

We want to thank to:

  Andreas Klöckner and everybody involved in the development of PyOpenCL 
  Ian Johnson (enjalot) for his excellent tutorial: "Adventures in PyOpenCL" 



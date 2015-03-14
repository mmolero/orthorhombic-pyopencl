if you use orthorhombic-pyopencl code in your research, we would appreciate the citation of the following article:

**"Accelerating numerical modeling of wave propagation through 2-D anisotropic materials using OpenCL",  Miguel Molero & Ursula Iturrarán-Viveros, Ultrasonics 53 (3), 2013, pages 815-822**

http://www.sciencedirect.com/science/article/pii/S0041624X12002612#



---


_Abstract_:

We present an implementation of the numerical modeling of elastic waves propagation, in 2D anisotropic materials, using the new parallel computing devices (PCD). Our study is aimed both to model laboratory experiments and explore the capabilities of the emerging PCDs by discussing performance issues.

> In the experiments a sample plate of an anisotropic material placed inside a water tank is rotated and, for every angle of rotation it is subjected to an ultrasonic wave (produced by a large source transducer) that propagates in the water and through the material producing some reflection and transmission signals that are recorded by a “point-like” receiver. This experiment is numerically modeled by running a finite difference code covering a set of angles θ ∈ [−50º , 50º], and recorded the signals for the transmission and reflection results. Transversely anisotropic and weakly orthorhombic materials are considered.

We accelerated the computation using an open-source toolkit called PyOpenCL, which lets one to easily access the OpenCL parallel computation API’s from high-level programming environment of Python. A speedup factor over 19 using the GPU is obtained when compared with the execution of the same program in parallel using a CPU multi-core (in this case we use the 4-cores that has the CPU). The performance for different graphic cards and operating systems is included together with the full 2-D finite difference code with PyOpenCL.


---


**Wave propagation for a weakly orthorhombic plate (rotated by an angle of 45º) in 2-D immersed in water tank.**

<a href='http://www.youtube.com/watch?feature=player_embedded&v=pw411BySLEY' target='_blank'><img src='http://img.youtube.com/vi/pw411BySLEY/0.jpg' width='425' height=344 /></a>





---


We want to thank to:

> Andreas Klöckner and everybody involved in the development of PyOpenCL

> Ian Johnson (enjalot) for his excellent tutorial: "Adventures in PyOpenCL"
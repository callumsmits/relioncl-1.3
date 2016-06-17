#RelionCL - RELION accelerated with GPUs using OpenCL
This is a port of Relion to accelerate key parts of the maximisation function using OpenCL. The code has been tested on a variety of Nvidia and AMD graphics cards successfully including:
- NVidia K80
- NVidia K20
- NVidia m2070
- AMD D700
- NVidia GTX 1080
- NVidia 780m

The software will detect and use all gpus within a MPI run and can operate with a heterogeneous mix of GPU's and instances running on the CPU. This allows the use of all available resources within a node. 

For optimal performance, run 2x more MPI instances on each compute node than GPUs on the node. Because not all code runs on the GPU, this ensures that the GPU is fully utilised.

All calculations are performed in double precision. If you wish to check the accuracy of your setup, uncomment (i.e. delete the //) line 7189 of ml_optimiser.cpp and run a refinement. The results of all calculations performed on the GPU will be compared to results from the CPU and any deviations will be reported. Additionally, by uncommenting line 7192 the speed of all calculations performed on the GPU and their relative speedup will be reported.

Please see the [RELION website](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page) for details on how to install and use RELION.

###FAQ
**Which GPU is best?**
Memory bandwith seems extremely important - more so than raw double precision throughput. For example, a GTX 1080 is roughly twice as fast as a K20 on the same data.

**What about other accelerators that can run OpenCL code?**
With minor modifications the code can run on other devices, but has been optimised for GPU use. It can run on the CPU, but this is pointless as the OpenCL code is less efficient than the original RELION code. It has been briefly tested using an Intel Xeon Phi and was able to run with minor modifications, but performance was poor. No testing has been done with an FPGA to date.

#RelionCL - RELION accelerated using OpenCL
This is a port of Relion to accelerate key parts of the maximisation function using OpenCL. The code has been tested on a variety of Nvidia and AMD graphics cards successfully including:
- NVidia K80
- NVidia K20
- NVidia m2070
- AMD D700
- NVidia GTX 1080
- NVidia 780m

The software will detect and use all gpus within a MPI run and can operate with a heterogeneous mix of GPU's and instances running on the CPU. This allows the use of all available resources. 

For optimal performance, run 2x more MPI instances on each compute node than GPUs on the node. Because not all code runs on the GPU, this ensures that the GPU is fully utilised.

Please see the [RELION website](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page) for details on how to install and use RELION.


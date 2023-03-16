# Meso-Micro-Simulation

#### Description
coupling the models of WRF to OpenFOAM to run a meso-micro scale simulation, contains the adapter to the WRF code and  various utilities of OpenFOAM

More specifically, the package contains the adapter to the WRF and the revised adapter to OpenFOAM, using the term of the PreCICE, which is the "glue" code library to couple the meso-scale and micro-scale simulation. 
For one thing, the adapter to the WRF model is a series patches, mainly adding and revising the code stored in the wrf/share directory.
For another thing, the adapter to OpenFOAM, whcih currently only supports OpenFOAM-7 but should work with other versions, extends to cover more variables output from WRF. In addition, the package contains a OpenFOAM utility to run a three-dimensional interpolation using the output from WRF to initialize the OpenFOAM simulation, which is expected to significantly reduces the time required for the OpenFOAM simulation to be compatible with the WRF simulation.

In general, the user is expected to download the meteorological data and intilaize the WRF simulation, which is usually a nested run, first. Afterwards, the WRF simulation within the finest grid is output to (a) initialize the OpenFOAM simualtion and (b) to the PreCICE library to further communicate with OpenFOAM.
Currently, only the one-way nesting is supported, as only the simulation results from the WRF side is ported to the OpenFOAM domain. The two-way nesting is still under devleopment.

#### Software Architecture
This software contains two parts,
(a) WRF-Adapter, which includes three different patch files to revise the source codes of WRF. More specifically, the patches target (1) the share/ directory, (2) the Registry/ directory and (3) the configure script
(b) OpenFOAM-Code, which includes three components could be built upon the installation of the OpenFOAM-7. More specifically, the components are (1) a library to do the three-dimensional interpolation, (2) a utility/exectuable to actually do the interpolation and (3) the revised OpenFOAM adapter 

#### Installation
The installation is divided according to the revisions made to the meso-scale model and micro-scale model.

For the meso-scale model:
1.  place the configure.patch, registry.patch and share.patch in the WRF-Adapter/ directory into the root directory of WRF source code where configure script is found
2.  change directory into the root directory of WRF source code
3.  patch -p0 < share.patch
4.  patch -p0 < registry.patch
5.  patch -p0 < configure.patch
6.  install WRF as usual (indicated by https://www2.mmm.ucar.edu/wrf/OnLineTutorial/). In the process of configuration, the script will prompt to ask whether to enable the PreCICE coupling
7.  when wrf.exe, real.exe appears in the main sub-directory of the root directory of WRF after compilation, the adapter is sucessfully built

For the micro-scale model:
1.  load the environment for the OpenFOAM-7, in the case the command wmake appears
2.  change directory into OpenFOAM-Code/threeDimensionalInterpolation
3.  wmake libso
4.  change directory into OpenFOAM-Code/volInteprolationSetField
5.  wmake
6.  change directory into OpenFOAM-Code/openfoam-adapter
7.  run Allwmake script
8.  check teh compilation log for any errors, and if trilinearInterpolateSetFields appears as exectuables


#### Instructions
For one want to run the Meso and micro scale coupled simulation, he/she needs to be farmiliar with both the WRF model and OpenFOAM model
A series of tutorial cases are under development, and the document is also under preparation

#### Contribution

1.  Fork the repository
2.  Create Feat_xxx branch
3.  Commit your code
4.  Create Pull Request


#### Gitee Feature

1.  You can use Readme\_XXX.md to support different languages, such as Readme\_en.md, Readme\_zh.md
2.  Gitee blog [blog.gitee.com](https://blog.gitee.com)
3.  Explore open source project [https://gitee.com/explore](https://gitee.com/explore)
4.  The most valuable open source project [GVP](https://gitee.com/gvp)
5.  The manual of Gitee [https://gitee.com/help](https://gitee.com/help)
6.  The most popular members  [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)

# aliqueit
[![Build status](https://ci.appveyor.com/api/projects/status/blfr3d71p8br2njo/branch/master?svg=true)](https://ci.appveyor.com/project/ChristianBeer/aliqueit/branch/master)

A C++ program for computing and verifying aliquot sequences. Designed to be run autonomously.

The actual factoring is mostly done with external programs. Trial division and Pollard rho are done internally, P-1/P+1/ECM are done with GMP-ECM, QS is done with Msieve or YAFU, and GNFS is done with ggnfs. Aliqueit tries to minimise the time spent factoring by running different methods to different depths depending on the composite size, though you're better off just using the "-y" argument and letting YAFU handle all that nowadays.

The code can be compiled on Linux using gcc and Windows using Visual Studio. You'll need a GMP lib for your platform.

Author: Mikael Klasson
Orignal site: http://mklasson.com/aliquot.php

Currently maintained by: Christian Beer

## Building on Linux
* Install [GMPlib](https://gmplib.org/) either using your package manager or compile from source on your own. Best place symlinks to header and library files in `3rdParty/` to your preferred version of GMP.
* If you have cmake installed you can change into `build/` and run `cmake .. && make` and get the executable. If you don't have cmake you can also use `make` in `src/` directly (after modifying paths in `src/Makefile`).

## Building on Windows
This is the way automated builds are working right now but if you install [cmake](https://cmake.org/) this also works locally.
* Get a precompiled GMP library from https://github.com/ShiftMediaProject/gmp/releases and unpack into `3rdParty/`. Or build your own GMP/MPIR library and copy the artifacts to `3rdParty/`.
* Generate the solution and project files using cmake like this: `cd build; cmake .. -G "Visual Studio 14 2015 Win64"`. Other versions should work too.
* Open the `build/aliqueit.sln` solution file and build the main project.

## External factoring programs

* [Most everything, incl. x64 versions](http://gilchrist.ca/jeff/factoring/)
* [GMP-ECM](http://gforge.inria.fr/projects/ecm/)
* [Msieve](http://sourceforge.net/projects/msieve/)
* [YAFU](http://sourceforge.net/projects/yafu/)
* [GGNFS](http://sourceforge.net/projects/ggnfs/)

## Tutorials

* [Set up Aliqueit on an Ubuntu computer](http://www.starreloaders.com/edhall/AliWin/AliqueitLinstall.html)
* [Set up Aliqueit and a GUI on a Windows computer](http://www.starreloaders.com/edhall/AliWin/AliWin.html)

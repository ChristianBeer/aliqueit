# aliqueit
A C++ program for computing and verifying aliquot sequences. Designed to be run autonomously.

The actual factoring is mostly done with external programs. Trial division and Pollard rho are done internally, P-1/P+1/ECM are done with GMP-ECM, QS is done with Msieve or YAFU, and GNFS is done with ggnfs. Aliqueit tries to minimise the time spent factoring by running different methods to different depths depending on the composite size, though you're better off just using the "-y" argument and letting YAFU handle all that nowadays.

Fairly portable Visual Studio 2008 source code is included in the windows package. You'll need a GMP lib if you want to compile the source.

Author: Mikael Klasson
Orignal site: http://mklasson.com/aliquot.php

## External factoring programs

* [Most everything, incl. x64 versions](http://gilchrist.ca/jeff/factoring/)
* [GMP-ECM](http://gforge.inria.fr/projects/ecm/)
* [Msieve](http://sourceforge.net/projects/msieve/)
* [YAFU](http://sourceforge.net/projects/yafu/)
* [GGNFS](http://sourceforge.net/projects/ggnfs/)

## Tutorials

* [Set up Aliqueit on an Ubuntu computer](http://www.starreloaders.com/edhall/AliWin/AliqueitLinstall.html)
* [Set up Aliqueit and a GUI on a Windows computer](http://www.starreloaders.com/edhall/AliWin/AliWin.html)

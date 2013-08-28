Installation {#install}
============

\page install

[git]:http://git-scm.com
[cmake]:http://cmake.org
[gcc]:http://gcc.gnu.org
[Ubuntu]:http://www.ubuntu.com

[TOC]

Obtaining the source code {#install_obtain}
=========================

The source code can be downloaded using the version control software [git].
You can obtain the current version of the software by typing the command

	git clone -b release_1.0 git://last.hit.bme.hu/toolbox/nihu.git

into the command line.
The command will create the directory `nihu` containing all the source files required for further steps.

Alternatively, you can download the source code as a packed archive.

Prerequisites {#install_prereq}
=============

In order to complie NiHu, the following prerequisites are needed:

In order to use the Matlab interface and compile mex files you must have:

- Matlab installed
- mex available

The installation process {#install_process}
========================

NiHu is installed using the free cross-platform make tool [cmake] on all supported platforms.
In the following the steps of the installation of the prerequisites and the compilation of the source code are discussed.

Installation on Unix systems {#install_unix}
----------------------------

This section presents how to install the prerequisites and compile NiHu from source code on Unix systems.
The example commands are given for and tested on [Ubuntu] 12.04.

### Installing [gcc]

Since NiHu requires a compiler that supports some features of the C++11 standard, you must ensure that you have such compiler on your system.
It is advised that you use the GNU Compiler Collection [gcc], which supports the required features from version 4.7.

You can install `gcc-4.7` if you have administrative rights on your computer in the following steps.

1. 	Add the test toolchain repository for `apt`
	
		sudo add-apt-repository ppa:ubuntu-toolchain-r/test
	
2. 	Update the list of available packages
	
		sudo apt-get update
	
3. 	Install `gcc-4.7` and `g++-4.7`

		sudo apt-get install gcc-4.7 g++-4.7

4. 	Add the newly installed compiler as an alternative to the old one.
	The example presents 

		sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.6 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.6 
		sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.7 40 --slave /usr/bin/g++ g++ /usr/bin/g++-4.7 

5.	Configure the system to use the newly installed compiler

		sudo update-alternatives --config gcc

	And select the appropriate choice.

6.	Make sure that the correct version is selected by

		gcc -v

### Installing [cmake]


Installation on Windows systems {#install_win}
-------------------------------

This is how to.



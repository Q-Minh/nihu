Installation {#install}
============

\page install

[git]:http://git-scm.com
[cmake]:http://cmake.org
[gcc]:http://gcc.gnu.org
[Ubuntu]:http://www.ubuntu.com
[doxygen]:http://www.stack.nl/~dimitri/doxygen
[TDM-GCC]:http://tdm-gcc.tdragon.net

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

Alternatively, or if you do not have administrative rights you can install [gcc] from source.

### Installing [cmake]

If you have administrative rights, you can easily install the newest version of [cmake] by the command

	sudo apt-get install cmake

### Configuration

For the installation you must first make a configuration step.
This is done by creating a directory for the build.
The build directory must not be the same as the source directory, however it can be a new subdirectory inside the main source directory.
In order to configure the installation you must run `cmake` inside the build directory.
The command `cmake` takes the path of the source directory as an argument.
The following example demonstrates the configuration inside the directory `nihu/build_dir`.
(Note: the source directory is located in the package at `nihu/src`.)

	cd nihu
	mkdir build_dir
	cd build_dir
	cmake ../src

For an advanced configuration, you can specify various options for the `cmake` command from the command line.
See [installation options](#install_cmake_options) for further details.

### Compiling the sources

All sources are compiled by calling the `make` command in the build directory.

	make

If you have [doxygen] installed, you can also compile the documentation by the command

	make doc

Installation on Windows systems {#install_win}
-------------------------------

### Installing gcc on windows

It is recommended that you compile NiHu on a Windows operating system also by using gcc.
A Windows version of `gcc-4.7` is available through the [TDM-GCC] project, you can download and install the binaries at their [download](http://tdm-gcc.tdragon.net/download) site.

### Installing cmake on windows

[cmake] is also available for windows, the executable can be downloaded [here](http://www.cmake.org/files/v2.8/cmake-2.8.11.2-win32-x86.exe).

### Configuration

For the proper configuration you must specify the option `-G "MinGW Makefiles"` for the `cmake` command.

	cd nihu
	md build_dir
	cd build_dir
	cmake ..\src -G "MinGW Makefiles"

### Compiling the sources

All sources are compiled by calling the `make` command in the build directory.

	mingw32-make.exe


Installation options for cmake {#install_cmake_options}
==============================

Eigen options {#install_eigen_options}
-------------



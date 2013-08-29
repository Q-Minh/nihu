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

1.  Add the test toolchain repository for `apt`
	
		sudo add-apt-repository ppa:ubuntu-toolchain-r/test

2.  Update the list of available packages
	
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

All sources are compiled by calling the `MinGW` version of the make command in the build directory.

	mingw32-make.exe

If you have [doxygen] installed, the documentation can be compiled by the command

	mingw32-make.exe doc

	
Configuration options {#install_cmake_options}
=====================

Installation options {#install_install_options}
--------------------

- **NIHU_INSTALL_DIR** Specifies the directory where the compiled binaries and header files are installed. The default install directory is the same as the build directory.

Eigen options {#install_eigen_options}
-------------

- **NIHU_EIGEN_INSTALL** When set as non-zero, the installer will not look for an installed version of Eigen, but Eigen headers are installed as a part of NiHu. This is the default option on Windows operating systems.
- **NIHU_EIGEN_DISABLE_PATCH** Disables the patching of Eigen files. The patch contains the removal of unnecessary `typedef`s in order to avoid compiler warnings of `gcc-4.8`.
- **NIHU_EIGEN_AS_PROJECT** When set to non-zero, the installation process will install Eigen as a standalone cmake project. This option only works on Unix system, due to the limitations of Eigen's cmake project. When this option is set all Eigen sources and tests are compiled.
- **NIHU_EIGEN_HEADERS_DIR** Specifies the location where Eigen files are copied.
- **NIHU_EIGEN_INSTALL_DIR** Specifies the location where Eigen header files are installed. You may have to set this parameter when using the `NIHU_EIGEN_AS_PROJECT` option since you may not have access rights for the default system directory.

Matlab options {#install_matlab_options}
--------------

- **NIHU_MATLAB_ROOT** Specifies where the installation will look for Matlab and its include directories. When not specified, the install process will look for Matlab in common directories, but it is possible that your Matlab installation will not be found by this search. In this case you should specify the Matlab root directory manually.
- **NIHU_FORCE_MEX_COMPILER** When set to non-zero, the compilation of mex files is done using the mex compiler of Matlab. Since the mex compiler invokes the system compiler, you must ensure that your system compiler supports the C++11 features. (This issue is mostly relevant in case of windows systems.)

Testing options {#install_testing_options}
---------------

- **NIHU_ENABLE_TESTING** When set to non-zero building of all tests are enabled.
- **NIHU_DISABLE_TEST_BUILD** When set to non-zero tests are excluded from the `make all` command. In this case tests can be compiled using separate `make` commands. This option is relevant only is `NIHU_ENABLE_TESTING` is set to non-zero.
- **NIHU_ENABLE_RUN_MATLAB_TESTS** When set to non-zero tests run from Matlab are included into all tests. Otherwise Matlab tests have to be run separately. (Note: this option is not recommended, since not all Matlab versions support command line mode, which is required for the tests in order to function properly.)

Documentation options {#install_doc_options}
---------------------

- **NIHU_HTML_DOC_PATH** When specified documentation will be created at the given path. Otherwise, the default path is `build_dir/doc/html`.
- **NIHU_MATHJAX_DISABLE** When set to non-zero MathJax is disabled in the documentation and formulae are displayed as images.
- **NIHU_MATHJAX_PATH** Specifies the path of the MathJax installation in order to make mathematical formulae appear with proper typeset in the documentation pages generated by Doxygen. Alternatively, MathJax is used from the online content distribution network (CDN) as the default option.
- **NIHU_ENABLE_DOC_INSTALL** When set to non-zero, the installation step for the documentation is also performed. Hence, the documentation will be a part of the resulting installation.

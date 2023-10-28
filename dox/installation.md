Installation {#install}
============

\page install Installation

[boost]:https://www.boost.org/
[boost-bin]:https://sourceforge.net/projects/boost/files/boost-binaries/
[CMake]:http://cmake.org
[cmake-download]:https://cmake.org/download/
[Doxygen]:http://www.stack.nl/~dimitri/doxygen
[Eigen]:http://eigen.tuxfamily.org
[FFTW]:http://www.fftw.org
[gcc]:http://gcc.gnu.org
[git]:http://git-scm.com
[gtest]:https://github.com/google/googletest
[Matlab]:http://www.mathworks.com/products/matlab
[Mex]:http://www.mathworks.com/help/matlab/create-mex-files.html
[Mingw64]:https://sourceforge.net/projects/mingw-w64/files/latest/download
[Ubuntu]:http://www.ubuntu.com
[VS]:https://visualstudio.microsoft.com/vs/features/cplusplus/

[tutorials]:\ref tutorials

[TOC]

Obtaining the source code {#install_obtain}
=========================

The source code of %NiHu can be downloaded using the version control software [git].
You can obtain the latest stable version of the software by typing the command

	git clone -b release_2.0 git://last.hit.bme.hu/toolbox/nihu.git

into the command line.
There is a nightly build available for %NiHu, which contains regular updates, however testing of the most recent features may not be complete.
You can download the nightly version similarly by the command

	git clone -b nightly git://last.hit.bme.hu/toolbox/nihu.git

The command will create the directory `nihu` containing all the source files required for further steps.
The directory containing the source files is referred to as _source directory_ in the following.
Inside the source directory, the folder `src` contains the C++ and Matlab source codes, whereas other folders contain the documentation and tools required for building the project.

Alternatively, you can download the source code as a packed archive.

Prerequisites and third party software {#install_prereq}
======================================

Build system
------------

To use %NiHu's automated, platform independent build system is based on the cross-platform tool [CMake][CMake].

On Linux, you can install `CMake` using tha pacakage manager, e.g.
	
	sudo apt install cmake

On Windows, you can download the [CMake installers][cmake-download].
\note The minimum required `CMake` version is `3.15`.

Tools required for compiling the source code
--------------------------------------------

To complie %NiHu, the following prerequisites are needed:

- A c++ compiler and linker that supports the C++14 standard. You will find further information on the compiler selection in the [next section](#install_process).
\note `NiHu` builds were tested using the [Gnu compiler collection][gcc], the [clang](http://clang.llvm.org/) compiler, and Microsoft [Visual Studio][VS].

- %NiHu relies on the template matrix library [Eigen][Eigen]. If you do not have Eigen installed on your computer, the installation process will download and install the necessary header files for the compilation of NiHu.
\note The current version of %NiHu was tested using `Eigen 3.3.7`.

- The evaluation of some advanced mathematical functions is based on the portable C++ library collection [Boost][boost].

	On Linux, you can install `Boost` using the package manager, e.g.
		
		sudo apt install libbboost-all-dev
		
	On Windows, you can install `Boost` binaries using the [Boost binary downloads][boost-bin], or you can build and install `Boost` manually.
	\note The minimum required `Boost` version is `1.67.0`.

- Fast Fourier transforms are evaluated using the C subroutine library [FFTW3][FFTW].

	On Linux, it is recommended that you install the FFTW3 package using the package manager, e.g.
		
		sudo apt install fftw3
		
	On Windows systems, if you do not have FFTW3 installed, the installation process will download FFTW automatically and create the necessary libraries for linking `NiHu` executables.
	\note The current version of %NiHu was tested using `FFTW 3.3.5`.
	\note On windows systems, it is recommended to add the FFTW3 directory to your `PATH` environment variable for easier usage of the libraries.
	
Optional components
-------------------
	
- %NiHu's test codes rely on [Googletest][gtest]. The source files of `gtest` are automatically downloaded during the installation procedure.

	\note The current version of %NiHu uses `gtest 1.7.0`.

- %NiHu's documentation is generated using [Doxygen][doxygen].

	On Linux, you can install doxygen using the package manager, e.g.
		
		sudo apt install doxygen
		

- Formulas in the documentation are produced using MathJax.


Matlab interface
----------------
		
To use %NiHu's `Matlab` interface and compile [mex][mex] files the following prerequisites are also needed:

- [Matlab] must be installed. Matlab versions 7.x, 8.x, and 9.x are supported.
- As a part of your Matlab installation you should also have the mex header files required to build C / C++ programs callable from Matlab.


Installation steps {#install_process}
==================

It is worth mentioning that since %NiHu is a template library, you do not need to compile any sources to use %NiHu's C++ core.
You can simply include the header files found in the directory `nihu/src` in order to compile your own C++ codes using the features implemented in %NiHu (see [below](#install_compile_cpp)).
However, the installation process lets you to compile %NiHu's libraries, tutorials, tests and `mex` files for %NiHu's `Matlab` interface.
Furthermore, if you complete the installation, you can make sure that you have all necessary prerequisites that are needed to use %NiHu.
Therefore, it is highly recommended to complete the installation before using the NiHu toolbox.

The installation is done in three steps as usual.

1. **Configuration**
	As a first step a directory for the build should be created referred to as _build directory_ in the following.
	\note The build directory _must not be the same_ as the source directory, however it can be a new subdirectory inside the main source directory.

2. **Build**
	The build step is the process, in which the source codes are compiled and the binaries of the executables and libraries are generated.

3. **Install**
	In the install step the generated binaries and includable header files are copied to the installation destination, called _installation directory_ in the follwing. By default this is the same as the build directory. When the installation process is complete, the build and source files are no longer necessary, they can be deleted if you do not want to use them.

In the following examples the build and installation directories will be located at `nihu/build_dir` and `nihu/install_dir`, respectively.
The steps of the installation of the prerequisites and the compilation of the source code are discussed for [Unix](#install_unix), [Windows with MinGW](#install_win_mingw), and [Windows with Visual Studio](#install_win_msvc) configurations in the sequel.
Alternatively, [Qt creator](#install_qt_creator) can also be used for importing the %NiHu project.

For an advanced configuration, you can specify various options for the `cmake` command from the command line.
See [installation options](#install_cmake_options) for further details.
Note that depending on your operating system and software configuration providing some command line options may be mandatory for a successful install.


Installation on Unix systems {#install_unix}
----------------------------

This section presents how to install the prerequisites and compile %NiHu from source code on Unix systems.

### Installing GCC

Since %NiHu requires a compiler that supports some features of the C++14 standard, you must ensure that you have such compiler on your system.
It is advised that you use the [GNU Compiler Collection][gcc], which supports the required features from version 6.1.

You can install the required build tools using the package manager, e.g.

	sudo apt install build-essentials

Make sure that the installed compiler version is at least 6.1

	gcc --version
	g++ --version

Alternatively, or if you do not have administrative rights, you can install [gcc] from source.

### Configuration

The configuration step is performed by calling the command `cmake` inside the build directory.
The command takes one command line argument, with the path to the directory of the C++ source files, which is `nihu/src` in our case.
The installation directory is defined by the option `-DNIHU_INSTALL_DIR="/path/to/install/dir"`.
The configuration is performed using the following commands.

	cd nihu
	mkdir build_dir
	cd build_dir
	cmake ../src -DNIHU_INSTALL_DIR="../install_dir"

### Compiling the sources

All sources are compiled by calling the `make` command in the build directory.

	make

If you have [doxygen] installed, you can also compile the documentation by the command

	make doc

### Installation steps

Finally, the binaries are installed using the command

	make install

With this step, the installation is completed and you can start using NiHu, see the [getting started](#install_get_started) section to see how to get started.

Installation on Windows systems using MinGW {#install_win_mingw}
-------------------------------------------

### Installing GCC on windows

You can compile %NiHu on a Windows operating system also by using gcc.
You can download a Windows version of `gcc-8` environment [here][Mingw64].

### Configuration

The configuration step is performed calling the command `cmake` inside the build directory.
The command takes one command line argument, with the path to the directory of the C++ source files, which is `nihu/src` in our case.
For the proper configuration you must specify the option `-G "MinGW Makefiles"` for the `cmake` command.
The installation directory is defined by the option `-DNIHU_INSTALL_DIR="/path/to/install_dir"`.

\note When specifying paths it is recommended to use the `/` (slash) character as a separator between subdirectories also on Windows systems. Usage of the `\` (backslash) can lead to misinterpretation of absolute and relative paths in cmake.

The configuration is performed by the following commands.
It is recommended that you use the `MinGW Command line` application for executing the following commands.

	cd nihu
	md build_dir
	cd build_dir
	cmake ../src -G "MinGW Makefiles" -DNIHU_INSTALL_DIR="../install_dir"
	
\note To be able to execute `cmake` from the command prompt, make sure that the `cmake` executable is included in your `PATH` environment variable.

### Compiling the sources

All sources are compiled by calling the `MinGW` version of the make command in the build directory.

	mingw32-make.exe

\note On 64-bit systems you can use the same `mingw32-make.exe` command to build 64-bit targets.

If you have [doxygen] installed, the documentation can be compiled by the command

	mingw32-make.exe doc

### Installation

Finally, the binaries are installed using the command

	mingw32-make.exe install

With this step, the installation is completed and you can start using NiHu, see the [getting started](#install_get_started) section to see how to get started.


Installation on Windows systems using Visual Studio {#install_win_msvc}
---------------------------------------------------

Later versions of [Visual Studio][VS] can also be used to compile %NiHu sources.
The following paragraphs demonstrate the process.

\note The current version of NiHu was tested using Visual Studio 2017 and 2019.

### Configuration

The configuration step is performed calling the command `cmake` inside the build directory.
The below example demonstrates how to set up the generator option of `cmake` for using Visual Studio.

	cd nihu 
	md build_dir
	cd build_dir
	cmake ..\src -G "Visual Studio 15 2017 Win64" -DNIHU_INSTALL_DIR="..\install_dir"
	
When the configuration is completed, various Visual Studio project and solution files are generated inside the build directory.
These projects can readily be opened 
	
### Build

To build %NiHu executables using Visual Studio, execute the following command in the command line.

	msbuild NiHu.sln /property:Configuration=Release

It is recommended to use the Developer Command Prompt for Visual Studio for executing the `msbuild` command.


Building the sources using Qt Creator {#install_qt_creator}
-------------------------------------

The %NiHu project can also be imported into Qt Creator by choosing `File` -> `Open project` and browsing the file `CMakeLists.txt` from the directory `nihu/src`.

Then, the project is built using the build tools configured inside Qt Creator.

	
Configuration options {#install_cmake_options}
=====================

You can define additional configuration options by using the `-D` option for the `cmake` command.
By these options you can specify parameters of `cmake` and NiHu itself.
All NiHu-specific options start with the prefix `NIHU_`.
The extra options are defined using the pattern

	cmake ../src -DNIHU_EXAMPLE_VARIABLE=1
	cmake ../src -DNIHU_EXAMPLE_PATH="path/to/example"

If you define more options at the same time, the name=value pairs should be separated by spaces.
You can find additional cmake-related configuration parameters on the website of [cmake].
Please note that it is not recommended to change the configuration settings that are not listed herein.

Compiler options {#install_compiler_options}
----------------

- **CMAKE_CXX_COMPILER** Specifies the compiler for C++ files for executing `make` commands. By default, the system default compiler is used. The value of this parameter should be the absolute full path to the executable file of the compiler. (For example, usage of the compiler g++-7 on Unix systems is achieved by the setting: `-DCMAKE_CXX_COMPILER="/usr/bin/g++-7"`)

FMM module options {#install_fmm_options}
------------------

The following options control the behavior of the Fast Multipole Method (FMM) module.

- **NIHU_DISABLE_FMM** Disables the FMM module, the FMM library and the FMM targets won't be built.
- **NIHU_FMM_DISABLE_PARALLEL** Set to non-zero to disable parallelization in fast multipole codes. This option is useful for writing manual parallel executions, such as freqeuncy parallelization, which may be more effective than parallelization for a single load case.
- **NIHU_FMM_TRAVERSE** This option controls the cluster tree traversal method and can be set either to `"BFS"` (Breadth-first search) or to `"DFS"` (Depth-first search).

Installation options {#install_install_options}
--------------------

- **NIHU_INSTALL_DIR** Specifies the directory where the compiled binaries and header files are installed. The default install directory is the same as the build directory. You can specify the installation directory either as a relative or a full path. You should verify that you have write access to the specified folder.

Eigen options {#install_eigen_options}
-------------

The following options control the setup of the matrix library Eigen during the installation process.
On Unix systems, the build system will automatically search for an existing installation of Eigen on your computer.
If the installed Eigen headers are found, these header files will be used for compiling %NiHu sources.
On Windows operating systems, %NiHu will not search for an existing Eigen installation, but Eigen headers are installed as a part of %NiHu.
You can override the default behavior by the parameters listed below.

- **NIHU_EIGEN_PATH** Specifies the full path to your existing Eigen installation, i.e., the path containing the directory `Eigen`. When this path is set, %NiHu will not search for an existing Eigen installation, but tries to use Eigen header files specified by this path.
- **NIHU_EIGEN_INSTALL** When set to a non-zero value, the installer will not look for an installed version of Eigen, but Eigen headers are installed as a part of NiHu. This is the default option on Windows operating systems. This option only has an effect when the path `NIHU_EIGEN_PATH` is not specified.
- **NIHU_EIGEN_ARCHIVE** Specifies the full path to a downloaded `Eigen` archive. NiHu will extract this archive instead of downloading the source from the predefined URL. This option only makes sense if the `NIHU_EIGEN_INSTALL` is set to non-zero.
\note If you experience a `cmake` or other error when downloading from the official URL, you should try this option.
- **NIHU_EIGEN_VERSION** When set and `Eigen` is installed as a part of NiHu (see above), you can specify which version to install. The default installation version is `3.3.7`. Please note that the currently supported versions are `3.2.0`, `3.2.7`, and `3.3.7`.

\note If you want to use Eigen as a stand-alone project, you can install it using cmake. Please consult the [Eigen website](http://eigen.tuxfamily.org) for further information.

FFTW options {#install_fftw_options} 
------------

The following options control the setup of the FFT library FFTW3.
On Unix systems, %NiHu will automatically search for an existing installation of FFTW3 on your computer.
On Windows operating systems, FFTW3 will be automatically downloaded and used as a part of %NiHu, furthermore, the FFTW static library files necessary for building %NiHu executables will have their own build target and will be built during the build process.
Alternatively, the variable `NIHU_FFTW_PATH` can be specified, and in this case the build system will look for the appropriate FFTW files in the specified folder.
You can override the default behavior using the parameters listed below:

- **NIHU_FFTW_PATH** Specifies the full path to an existing FFTW directory. The directory should contain the fftw header file `fftw3.h` and the libraries.
- **NIHU_FFTW_INSTALL** When set to a non-zero value, the installer will not look for an installed version of FFTW, but FFTW headers and libraries are installed as a part of NiHu. This is the default option on Windows operating systems, unless the variable `NIHU_FFTW_PATH` is specified.
- **NIHU_FFTW_ARCHIVE** Specifies the path to a downloaded FFTW3 archive. The build system will extract this archive instead of downloading the source from the predefined URL. This option only makes sense if the variable `NIHU_FFTW_INSTALL` is set to non-zero.
- **NIHU_FFTW_VERSION** Specifies which FFTW version to download. Currently only version `3.3.5` is supported.

Boost options {#install_boost_options}
-------------

During the configuration step the build system will automatically look for Boost headers and libraries.
If your Boost installation is not found automatically, try to specify the option `NIHU_BOOST_PATH` as explained below.
The default behavior of finding the Boost libraries can be overridden using the following parameters:

- **NIHU_BOOST_PATH** Specifies the full path to the boost installation, should be the path where the directory named `boost` is found.
- **NIHU_BOOST_LIB_PATH** Specifies the full path to the boost compiled libraries (necessary on windows)

Matlab options {#install_matlab_options}
--------------

In order to use the Matlab interface of %NiHu and to compile `mex` source files the setup process must find an existing Matlab installation on your computer.
%NiHu will search for the root directory of your Matlab installation and also for the `mex` C++ header files.
You can customise the Matlab related settings of the installation process by the parameters listed below.

- **NIHU_MATLAB_PATH** Specifies where the installation will look for Matlab and its include directories. When not specified, the install process will look for Matlab in common directories, but it is possible that your Matlab installation will not be found by this search. In this case you should specify the Matlab path manually.
- **NIHU_MATLAB_FORCE_MEX_COMPILER** When set to non-zero, the compilation of mex files is done using the mex compiler of Matlab. Since the mex compiler invokes the system compiler, you must ensure that your system compiler supports the necessary features of the C++14 standard.

Testing options {#install_testing_options}
---------------

NiHu comes with various test sources, such as
- unit tests for testing the functionalities of each module separately
- numerical tests for validating the computations on academic examples
- Matlab tests for testing mex interface

By default, the tests are excluded from the build process, however, you can turn on testing and control the build and installation parameters of all tests using the following settings.

- **NIHU_ENABLE_TESTING** When set to a non-zero value building of all tests are enabled.
- **NIHU_ENABLE_TEST_INSTALL** When set to non-zero the test executables are included in the installation. This option is only relevant if `NIHU_ENABLE_TESTING` is set to non-zero.
- **NIHU_DISABLE_TEST_BUILD** When set to non-zero tests are excluded from the `make all` command. In this case tests can be compiled using separate `make` commands. This option is only relevant if `NIHU_ENABLE_TESTING` is set to non-zero.
- **NIHU_ENABLE_RUN_MATLAB_TESTS** When set to non-zero tests run from Matlab are included into all tests. Otherwise Matlab tests have to be run separately by using the test script provided with the compiled mex files. This option is only relevant if `NIHU_ENABLE_TESTING` is set to non-zero.
\note This option is not recommended, since not all Matlab versions support command line mode, which is required for the tests in order to function properly.)

Tutorial options {#install_tutorial_options}
----------------

You can control the compilation of tutorials by the following parameters.

- **NIHU_DISABLE_TUTORIALS** When set to a non-zero value building of all tutorials are disabled.


Documentation options {#install_doc_options}
---------------------

If you have Doxygen installed, you can build the documentation of NiHu on your own computer and access all help and tutorial documentation locally.
The installation process will automatically look for an existing installation of Doxygen on your computer.
It is recommended to you use [Doxygen] version up from 1.8.4. for compiling the documentation.

- **NIHU_HTML_DOC_DIR** When specified documentation will be created at the given path. Otherwise, the default path is `build_dir/doc/html`.
- **NIHU_MATHJAX_DISABLE** When set to non-zero MathJax is disabled in the documentation and formulae are displayed as images.
- **NIHU_MATHJAX_PATH** Specifies the path of the MathJax installation in order to make mathematical formulae appear with proper typeset in the documentation pages generated by Doxygen. Alternatively, MathJax is used from the online content distribution network (CDN) as the default option. This option is only relevant when `NIHU_MATHJAX_DISABLE` is not specified or set to zero.
- **NIHU_ENABLE_DOC_INSTALL** When set to non-zero, the installation step for the documentation is also performed. Hence, the documentation will be a part of the resulting installation. The html documentation is installed in the folder `install_dir/doc`.

\note If you want to rebuild and reinstall NiHu with changing the configuration options, it is advised that you clean the build and install directories before doing so. Not cleaning these directories may result in unwanted options stored in the `cmake` cache.

Getting started {#install_get_started}
===============

After the installation process is successfully completed you should be ready to use NiHu on your computer. The brief description in the following will help you with [running NiHu tests](#install_run_tests), [compiling your own C++ sources](#install_compile_cpp) and [starting to use the Matlab interface](#install_matlab_interface).

Running tests {#install_run_tests}
-------------

If you have built the NiHu test executables (by setting the option `-DNIHU_ENABLE_TESTING=1` for the `cmake` command), you can run them by executing the command `ctest` inside the build directory. (For other testing related options please refer to [Testing options](#install_testing_options).)

	cd build_dir
	ctest
	
The command `ctest` automatically executes all tests and validates the output results. You should see that all tests are passed. If you want to execute only a part of the generated tests, you can execute the command `ctest` in any of the subdirectories of the folder `test` inside your build directory. You can also execute the tests one by one in order to get the details of the test outputs.

Compiling C++ sources using NiHu {#install_compile_cpp}
--------------------------------

In order to create and compile your own C++ sources relying on the toolbox, you only have to include the installed header files from the installation directory, `install_dir/include`.
It is convienient to include the NiHu headers in the following manner.

~~~~~~~~{.cpp}
#include "nihu/core/weighted_residual.hpp"
#include "nihu/tmp/algorithm.hpp"
...

int main(void)
{
	...
	return 0;
}
~~~~~~~~

If your source file is named `example.cpp` you can compile it and create the executable `example` using the `g++` compiler with specifying the include directories by the `-I` switch using the command

	g++ example.cpp -std=c++14 -I/path/to/nihu_install_dir/include -o example
	
By using the above pattern you only have to add one directory to the include path definitions of the compiler, as demonstrated above.
\note You should always use the `-std=c++14` option when compiling C++ sources using %NiHu.
	
Using the Matlab interface {#install_matlab_interface}
--------------------------

In order to use the functions of the Matlab interface, you should add the Matlab interface path of your %NiHu installation to the search path of Matlab. This can be done by executing the `install` command inside the `matlab` folder of your NiHu installation.

	>> cd '/path/to/nihu_install/matlab'
	>> install

This will add the necessary directories to Matlab's search path and enables you to call %NiHu's Matlab and mex functions from any of your Matlab scripts. After the `install` script is successfully executed you will see a window displaying the demo applications for the Matlab interface of NiHu. You can open these demos anytime by the command

	>> demo toolbox nihu

\note There is an issue that can occur when executing mex files on a Unix systems. Matlab can give you an error that says `libstdc++.so.6: version GLIBCXX_...` is not found when executing a mex file generated by the %NiHu installation process. This is because Matlab not necessarily support newer versions of `gcc`, see [supported compilers](https://www.mathworks.com/support/requirements/supported-compilers.html?sec=glnxa64). You can fix this issue by redirecting the soft link `libstdc++.so.6` in the directory `MATLAB_ROOT/sys/os/glnxa64/` to the current version of the stdc++ library. This can be achieved by executing the following commands

		cd /path/to/MATLAB/sys/os/glnxa64
		sudo mv libstdc++.so.6 libstdc++.so.6.old
		sudo ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.17 libstdc++.so.6

Replacing the soft link should fix the problem and `mex` files should work from this point after restarting Matlab.
You can find a longer discussion of this problem on [stackoverflow](http://stackoverflow.com/questions/17000903/mex-compiling-on-64-bit-linux-usr-bin-ld-cannot-find-lstdc).

\note On Windows, with MinGW build, older Matlab versions may report the error "Invalid MEX-file <i>file name</i>: The specified procedure could not be found."
If this happens, copy all the .dll files of your MinGW installation (located at MinGW/mingw64/bin) to MATLAB/bin/win64 directory.
Note that it is useful to create a backup of libstdc++-6.dll first.

Further steps {#install_further_steps}
-------------

Following this documentation you should be able to compile and install %NiHu from source.
You should also be able to compile your own C++ files using the header files of %NiHu or create and run own Matlab files using %NiHu's Matlab interface.

In order to find introductory examples of applying NiHu for basic engineering problems, please see the [tutorials] that will guide you through further steps.

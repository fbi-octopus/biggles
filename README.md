Getting Started
===============

The Biggles tracker is written in C++ under Ubuntu 16.04. and requires the following libraries be installed on your system:

-   [Boost](http://boost.org/)
-   [Eigen3](http://eigen.tuxfamily.org/)
-   [Jansson](http://www.digip.org/jansson/)

In addition you will need the following tools/utility programs:

-   [CMake \>= 2.8](http://cmake.org/)
-   Boost.Python and the Python development libraries (optional)
-   LaTeX (optional)
-   [Doxygen](http://doxygen.org/) (optional)
-   [Sphinx](http://sphinx-doc.org/) (optional)

Building
--------

It is recommended that you build the source out of tree. Checkout the Biggles tracker from the source control, make a build directory and run CMake:

```console
$ git clone git@github.com:fbi-octopus/biggles
$ cd biggles
$ mkdir build
$ pushd build; cmake -DCMAKE_INSTALL_PREFIX=~/.local ..; popd
$ make install test
$ make doc # If you have Doxygen and/or sphinx installed
```

If there are any libraries missing the `cmake` call above will report which libraries need to be installed.

All tests *should* pass. Replace \~/.local with the directory you want to install Biggles to. If it is \~/.local, make sure the \~/.local/bin is on your path and \~/.local/lib is in your library path. To build a project using Biggles, you may need to set the PKG\_CONFIG\_PATH environment variable to \~/.local/lib/pkgconfig.

Documentation
-------------

If Doxygen was installed then there should be some generated C++ documentation in `build/doc/html`. If sphinx was installed and the Python bindings are built then there is documentation about them in `build/python/doc/html`.

Tracking data
-------------

Example data can be found [here](http://dx.doi.org/10.5286/edata/733). The basic command is: 

```console
$ biggles track --dual input.json
```

Additional options can be displayed by

```console
$ biggles track --help
```


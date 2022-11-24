This is the first submission of this package.

## Test environments

* local Windows 10 install, R 4.1
* GitHub via Actions:
  * Mac OS X 10.15.7, R 4.1
  * Ubuntu 20.04.2 LTS, R 4.1
  * Ubuntu 20.04.2 LTS, R devel
  * Microsoft Windows Server 2019 10.0.17763 Datacenter
* R-Hub
  * Windows Server 2008 R2 SP1, R 4.1, 32/64 bit
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Debian Linux, R-devel, GCC ASAN/UBSAN
* win-builder
  * R version 4.1.0 (2021-05-18) using platform x86_64-w64-mingw32 (64-bit)

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs on Mac OS and Ubuntu.

### On Windows there is one WARNING due to compiler warnings.

These warnings come from the Boost headers. The Boost maintainers have evaluated them, consider them false positives, and used #pragmas to disable the warnings; see here:
https://github.com/boostorg/container/commit/6504af87080ec0f5193e0cd623795dedc4a5d9c3

The maintainer of the BH R package re-enabled the compiler warnings by commenting out the #pragmas here:
https://github.com/eddelbuettel/bh/commit/54182166369ef0ac1e7a58ef331afc02f1c5dd2c

As a result, my code generates these false-positive warnings:

* checking whether package 'rtree' can be installed ... WARNING

Found the following significant warnings:

  D:/a/_temp/Library/BH/include/boost/container/detail/copy_move_algo.hpp:183:19: warning: 'void* memmove(void*, const void*, size_t)' writing to an object of type 'value_type' {aka 'struct std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, __gnu_cxx::__normal_iterator<const std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, unsigned int>*, std::vector<std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, unsigned int> > > >'} with no trivial copy-assignment; use copy-assignment or copy-initialization instead [-Wclass-memaccess]

  D:/a/_temp/Library/BH/include/boost/container/detail/copy_move_algo.hpp:212:19: warning: 'void* memmove(void*, const void*, size_t)' writing to an object of type 'boost::movelib::detail::iterator_to_element_ptr<std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, __gnu_cxx::__normal_iterator<const std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, unsigned int>*, std::vector<std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, unsigned int> > > >*>::element_type' {aka 'struct std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, __gnu_cxx::__normal_iterator<const std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, unsigned int>*, std::vector<std::pair<boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>, unsigned int> > > >'} with no trivial copy-assignment; use copy-assignment or copy-initialization instead [-Wclass-memaccess]

### On local Windows there is one additional NOTE

This also appears to be a false positive. There are no calls to 'abort', 'exit' or 'printf' in my code.

> checking compiled code ... NOTE
  Note: information on .o files for x64 is not available
  File 'C:/Research/rtree.Rcheck/rtree/libs/x64/rtree.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  
  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.


## Downstream dependencies

This is a first release to CRAN, there are currently no downstream dependencies for this package.

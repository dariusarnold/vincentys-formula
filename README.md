# vincentys-formula
Speed comparison between Python (Numpy and pure Python), C and C++.

This repository contains Vincenty's formula in
  - C, to be used with the ctypes Python module.
  - C++, compiled with pybind11
  - Python (Numpy and pure Python)
  
The formula calculates the azimuth and distance between two points on the surface of a spheroid.
For more information about the formula see: https://en.wikipedia.org/wiki/Vincenty's_formulae

## License and usage

Licensed under GPLv3.
If you use the program for something cool, I'd love to hear about it!
Please report any problems, ideas, etc. as an issue.

### Compile C code to a shared library

CMake is the preferred way to build and install the library:
```
cmake -B build -S .
cmake --build build --target vinc
cmake --install build --prefix ./install
```

Alternatively you could build directly with GCC:
```
gcc -O3 -shared -lm -Wl,-soname,vinc -o vinc.so -fPIC -Isrc/ext src/vinc.c
```

### Compile c++ code to an extension module

The vinc_cpp Python extension module can be built by CMake like the vinc C library. 

If you insist, you can build it manually with (adapt to local paths or versions):
```
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) src/vinc.cpp -Isrc/ext -o vinc_cpp$(python3-config --extension-suffix) -I /usr/include/python3.11
```

### Generate test results

The script requires you to build vinc and vinc_cpp. 
Execute the vinc_c_vs_python.py file, see --help for options:
```
./vinc_c_vs_python.py --help
usage: VincentyCompare [-h] VINC_SO VINC_CPP

Vincenty Comparison Python vs C

positional arguments:
  VINC_SO     Path to vinc shared library file to use, typically libvinc.so.
  VINC_CPP    Path to C++ Python extension module, for example vinc_cpp.cpython-311-x86_64-linux-gnu.so.

options:
  -h, --help  show this help message and exit
```

On a Xeon E3-1246v3 @ 3.50GHz
```
Testing function vinc
NumPy best average out of 1000 runs:  9.9E-05 secs
Python best average out of 1000 runs: 9.7E-06 secs
C best average out of 1000 runs:      1.1E-06 secs
C++ best average out of 1000 runs:    5.8E-07 secs

Numpy results:  559580.5062148898 -3.0210528410061204
Python results: 559580.5062148898 -3.0210528410061204
C results:      559567.3915137552 -3.021053212988765
C++ results:    559580.5061678728 -3.0210528410061204
```

## Conclusion

Numpy is slower than pure Python since only scalars are used during the calculation.
C takes 1.1 µs per function call, C++ 0.6 µs and python 10 µs. Binding the extension with pybind11 requires less boilerplate code than ctypes, and is not less efficient. The generated module is larger though (122 kB in C++ vs 16 kB for C).

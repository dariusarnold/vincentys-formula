# vincentys-formula
Vincenty's formula for distance in C, to be used with the ctypes Python module.
The formula calculates the azimuth and distance between two points on the surface of a spheroid.
For more information about the formula see: https://en.wikipedia.org/wiki/Vincenty's_formulae

Vincenty's formula is also implemented in Python with NumPy and pure Python math functions for a speed comparison.

## License and usage

Licensed under GPLv3.
If you use the program for something cool, I'd love to hear about it!
Please report any problems, ideas, etc. as an issue.

### Compile to a shared library 

CMake is the preferred way to build the library:
```
cmake -B build -S .
cmake --build build --target vinc
```

Alternatively you could build directly with GCC:
```
gcc -shared -lm -Wl,-soname,vinc -o vinc.so -fPIC -Isrc/ext src/vinc.c
```

### How to use in Python
See the vinc_c_vs_python.py file help:

```
python3 vinc_c_vs_python.py --help
usage: VincentyCompare [-h] VINC_SO

Vincenty Comparison Python vs C

positional arguments:
  VINC_SO     Path to vinc shared library file to use.

options:
  -h, --help  show this help message and exit
```

### Generate test results
Execute the vinc_c_vs_python.py file

```
> python3 vinc_c_vs_python.py build/libvinc.so
Trying to load /home/darius/git/vincentys-formula/cmake-build-release/libvinc.so

Testing function vinc with Numpy Math
Python best average out of 1000 runs: 5.861742600246E-05 secs
C best average out of 1000 runs:      8.906520015444E-07 secs
C is 65.81 times faster
Python results: 559580.5062148898 -3.0210528410061204
C results:      559580.5061678728 -3.0210528410061204
Difference distance: 8.40219627E-09%
Difference azimut:   0.00000000E+00%

Testing function vinc with Python Math
Python best average out of 1000 runs: 5.472176999319E-06 secs
C best average out of 1000 runs:      8.970049966592E-07 secs
C is 6.10 times faster
Python results: 559580.5062148898 -3.0210528410061204
C results:      559580.5061678728 -3.0210528410061204
Difference distance: 8.40219627E-09%
Difference azimut:   0.00000000E+00%

Testing function trans
Python best average out of 1000 runs: 6.219373200292E-05 secs
C best average out of 1000 runs:      9.205289970851E-07 secs
C is 67.56 times faster
Python results: -554806.1287020878 -67202.02083485905
C results:      -554806.1286555919 -67202.02082922714
Difference distance: 8.38056735E-09%
Difference azimut:   8.38056735E-09%
```

# vincentys-formula
Vincenty's formula for distance in C, to be used with the ctypes Python module.

Also implemented in Python with NumPy and pure Python math functions for a speed comparison.

### Compile to a shared library 
gcc -shared -lm -Wl,-soname,vinc -o vinc.so -fPIC main.c

### How to use in Python
See the vinc_c_vs_python.py file

### Generate test results
Execute the vinc_c_vs_python.py file

```
Testing function vinc with Numpy Math
Python best average out of 1000 runs: 6.780600547791E-05 secs
C best average out of 1000 runs:      1.336097717285E-06 secs
C is 50.75 times faster
Python results: 559580.506168 -3.02105284101
C results:      559567.391514 -3.02105321299
Difference distance: 2.34371307E-03%
Difference azimut:   -1.23130120E-05%

Testing function vinc with Python Math
Python best average out of 1000 runs: 1.299405097961E-05 secs
C best average out of 1000 runs:      1.357078552246E-06 secs
C is 9.58 times faster
Python results: 559580.506168 -3.02105284101
C results:      559567.391514 -3.02105321299
Difference distance: 2.34371307E-03%
Difference azimut:   -1.23130120E-05%

Testing function trans
Python best average out of 1000 runs: 6.957101821899E-05 secs
C best average out of 1000 runs:      1.410961151123E-06 secs
C is 49.31 times faster
Python results: -554806.128656 -67202.0208292
C results:      -554793.184347 -67200.2435222
Difference distance: 2.33317734E-03%
Difference azimut:   2.64479254E-03%
```

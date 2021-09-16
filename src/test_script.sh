export LD_LIBRARY_PATH=/usr/local/gcc-6.3.0/lib64
make
gfortran-6 -o math_test.exe math_test.f90 math.o /usr/local/lapack-3.10.0-gcc6/libblas.a /usr/local/lapack-3.10.0-gcc6/liblapack.a
./math_test.exe

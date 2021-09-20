export LD_LIBRARY_PATH=/usr/local/gcc-6.3.0/lib64 /usr/local/lapack-3.10.0-gcc6/libblas.a /usr/local/lapack-3.10.0-gcc6/liblapack.a
make
gfortran-6 -o prim_test.exe prim_test.f90 math.o /usr/local/lapack-3.10.0-gcc6/libblas.a /usr/local/lapack-3.10.0-gcc6/liblapack.a -lblas -llapack
./math_test.exe

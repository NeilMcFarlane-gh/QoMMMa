export LD_LIBRARY_PATH=/../../../../../usr/local/gcc-6.3.0/lib64
make
gfortran-6 -o math_test.exe math_test.f90 math.o
./math_test.exe

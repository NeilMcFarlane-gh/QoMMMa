FF=gfortran-6
LD=gfortran-6
FFLAGS=
LIBS=
LFLAGS=
RM=rm -f
ECHO=echo

BIN=../../bin

FREQ=$(BIN)/freq.x
FREQ_OBJ= \
calculate_freqs.o Jacobi_diag.o Jacobi.o projection.o read_dimensions.o ReadInput.o \
write_hessian.o

FREQ_PRJ=$(BIN)/prj_freq.x
FREQ_PRJ_OBJ= \
calculate_prj_freqs.o Jacobi_diag.o Jacobi.o projection.o read_dimensions.o ReadInput_prj.o \
write_hessian.o  

all:	$(FREQ) $(FREQ_PRJ)

$(FREQ) : $(FREQ_OBJ)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

$(FREQ_PRJ) : $(FREQ_PRJ_OBJ)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

%.o : %.f90
	$(FF) $(FFLAGS) -c  $< -o $@ 

clean:
	$(ECHO) Cleaning up...
	$(RM) *.o *.mod *~


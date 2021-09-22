FF=gfortran
LD=gfortran
FFLAGS=
LIBS=
LFLAGS=
RM=rm 
ECHO=echo

all:
	make -C src FF=$(FF) LD=$(LD) FFLAGS=$(FFLAGS) LIBS=$(LIBS) LFLAGS=$(LFLAGS) RM=$(RM) ECHO=$(ECHO)
	make -C src/tinker FF=$(FF) LD=$(LD) FFLAGS=$(FFLAGS) LIBS=$(LIBS) LFLAGS=$(LFLAGS) RM=$(RM) ECHO=$(ECHO)
	make -C src/freq FF=$(FF) LD=$(LD) FFLAGS=$(FFLAGS) LIBS=$(LIBS) LFLAGS=$(LFLAGS) RM=$(RM) ECHO=$(ECHO)
	./qomdoc.py html
        
clean:
	make -C src clean FF=$(FF) LD=$(LD) FFLAGS=$(FFLAGS) LIBS=$(LIBS) LFLAGS=$(LFLAGS) RM=$(RM) ECHO=$(ECHO)
	make -C src/tinker clean FF=$(FF) LD=$(LD) FFLAGS=$(FFLAGS) LIBS=$(LIBS) LFLAGS=$(LFLAGS) RM=$(RM) ECHO=$(ECHO)
	make -C src/freq clean FF=$(FF) LD=$(LD) FFLAGS=$(FFLAGS) LIBS=$(LIBS) LFLAGS=$(LFLAGS) RM=$(RM) ECHO=$(ECHO)


FF=gfortran-6
LD=gfortran-6
FFLAGS=-fbounds-check #-Wall -Wextra -Wno-tabs
LIBS=/usr/local/lapack-3.10.0-gcc6/libblas.a /usr/local/lapack-3.10.0-gcc6/liblapack.a
LFLAGS=-llapack -lblas
RM=rm -f
ECHO=echo

BIN=../bin

QOMMMA=$(BIN)/qommma8.x
QOMMMA_OBJ= \
nrtype.o  coordinates.o  optimdata.o \
math.o primitive.o dlc_constraint.o delocalised.o \
mainprogram.o \
alloc_optim.o read_converg.o alloc_coord.o \
read_checkfile.o  read_expl_coord.o  read_tinkergrad.o \
read_ab_initio.o read_qmlatgrad.o read_mulliken.o \
evaluate_overall_grad.o update_opt_geometry.o update_conj_opt.o \
update_full_geometry.o evaluate_spring_force.o check_constrained_atoms.o \
evaluate_constraint_cart.o write_gradcorrection.o write_report.o write_checkfile.o \
evaluate_constraint_dlc.o det_sm_coord.o  write_sm_coord.o  write_hessgeom.o  write_mulliken.o \
disp_corr.o write_prim_coord.o

POSTMICRO=$(BIN)/postmicro8.x
POSTMICRO_OBJ= \
nrtype.o  coordinates.o  optimdata.o \
math.o primitive.o delocalised.o dlc_constraint.o  \
process_postmicro.o alloc_coord.o  alloc_optim.o read_expl_coord.o  read_hessgeom.o \
short_read_checkfile.o collategeoms.o write_coord.o write_charges.o write_resp.o \
check_constrained_atoms.o

SETUP=$(BIN)/setup8.x
SETUP_OBJ= \
nrtype.o  coordinates.o optimdata.o  \
math.o primitive.o dlc_constraint.o delocalised.o \
alloc_optim.o alloc_coord.o  read_expl_coord.o  write_first_report.o  read_charges.o  \
write_charges.o   write_resp.o   write_inactive.o  write_first_checkfile.o  \
read_initfile.o  setup_first.o det_sm_coord.o write_sm_coord.o write_atomic.o \
check_constrained_atoms.o initialize_update_geom.o initialize_prims.o \
initialize_step.o update_step_geometry.o

QOMMMA_MECP=$(BIN)/qommma_mecp.x
QOMMMA_MECP_OBJ= \
nrtype.o  coordinates.o  optimdata.o \
math.o primitive.o delocalised.o dlc_constraint.o  \
mainprogram_mecp.o \
alloc_optim_mecp.o  read_converg.o alloc_coord.o \
read_checkfile_mecp.o  read_expl_coord.o  read_tinkergrad.o \
read_ab_initio_mecp.o Effective_Gradient_mecp.o read_qmlatgrad_mecp.o read_mulliken_mecp.o \
evaluate_overall_grad_mecp.o update_opt_geometry_mecp.o update_conj_opt.o \
update_full_geometry.o evaluate_spring_force.o check_constrained_atoms.o \
evaluate_constraint_cart.o write_gradcorrection.o write_report_mecp.o write_checkfile_mecp.o \
evaluate_constraint_dlc.o det_sm_coord.o  write_sm_coord.o  write_hessgeom.o  write_mulliken.o \

GSM_REPARAM=$(BIN)/gsm_reparam.x
GSM_REPARAM_OBJ= \
nrtype.o  coordinates.o  optimdata.o \
math.o primitive.o delocalised.o dlc_constraint.o \
reparam_string.o \
alloc_coord.o alloc_optim.o read_expl_coord.o \
short_read_checkfile.o initialize_step.o update_step_geometry.o \
write_prim_coord.o check_constrained_atoms.o \



all:	$(QOMMMA) $(POSTMICRO) $(SETUP) $(GSM_REPARAM)

$(QOMMMA) : $(QOMMMA_OBJ)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

$(POSTMICRO) : $(POSTMICRO_OBJ)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

$(SETUP) : $(SETUP_OBJ)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

#$(QOMMMA_MECP) : $(QOMMMA_MECP_OBJ)
	#$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

$(GSM_REPARAM) : $(GSM_REPARAM_OBJ)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

%.o : %.f90
	$(FF) $(FFLAGS) -c  $< -o $@ 

clean:
	$(ECHO) Cleaning up...
	$(RM) *.o *.mod *~

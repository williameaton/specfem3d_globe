# includes default Makefile from previous configuration
include Makefile

# test target
default: test_sinsq_stf

## compilation directories
O := ./obj

OBJECTS = \
	$(specfem3D_SOLVER_OBJECTS) \
	$(specfem3D_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

test_sinsq_stf:
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_sinsq_stf test_sinsq_stf.f90 -I./obj $(OBJECTS) $(MPILIBS)


# includes default Makefile from previous configuration
include Makefile

# test target
default: test_siem_mesh

## compilation directories
O := ./obj

OBJECTS = \
	$O/meshfem3D_par.check_module.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

test_siem_mesh:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_siem_mesh test_siem_mesh.f90 -I./obj $(OBJECTS)


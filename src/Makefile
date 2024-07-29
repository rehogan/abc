#!/bin/sh  -V
# gfortran -Wimplicit-interface abc.f90 -o abc.exe
FC = gfortran
FFLAGS = -g -Wimplicit-interface
EXE = abc_db.exe
LIBS =
LDFLAGS = 

%.o: %.f90 
	@(echo $(FC) $(FFLAGS) -c $< )
	@$(FC) $(FFLAGS) -c $<
	
MODS =  mod_precision.f90 \
	mod_constants.f90 \
	mod_air_properties.f90 \
	mod_air_table.f90 \
	mod_bc.f90 \
	mod_files.f90 \
	mod_global_parameters.f90 \
	mod_interface_definitions.f90 \
	mod_lin_interp.f90 \
	mod_mat_prop.f90 \
	mod_nodes_elements.f90 \
	mod_parsing.f90 \
	mod_rad_encl.f90 \
	mod_rad_encl_lib.f90 \
	mod_time_params.f90 \
	mod_to_lower.f90

SRCS = abc.f90 \
	air_props.f90 \
	bc_del_natl_conv.f90 \
	bc_del_ff_rad.f90 \
	bc_del_nl_flux.f90 \
	bc_del_spec_flux.f90 \
	bc_del_st.f90 \
	band.f90 \
	compare_string.f90 \
	cond_cap.f90 \
	c_prop.f90 \
	el_matrix_to_global_matrix.f90  \
	el_bc_matrix_to_global_matrix.f90 \
	elem_cap_quad.f90 \
	elem_cond_quad_iso.f90 \
	find.f90 \
	find_string.f90 \
	flux_natl_conv.f90 \
	get_en_cnt.f90 \
	get_prop.f90 \
	initialize.f90 \
	mat_vec_mult.f90 \
	natl_conv_corr.f90 \
	open_files.f90 \
	read_aba_node_elem_blk.f90 \
	read_conv_bc.f90 \
	read_encl_rad_bc.f90 \
	read_glob_param.f90 \
	read_input.f90 \
	read_record.f90 \
	read_spec_t_bc.f90 \
	read_table.f90 \
	read_therm_prop.f90 \
	read_title.f90 \
	record_type.f90 \
	strip_blank.f90 \
	upcase.f90 \
	verify_lib.f90 \
	vol_of_scv_quad.f90

OBJS = $(MODS:f90=o) $(SRCS:f90=o)

# the default target is "abc_db.exe" which is in the ".depend" file
# So, make without a target, "make" will default to the target "abc_db.exe" 

# include the ".depend" file.
# It contains the list of objects, the executable target, and the dependencies 
include .depend

depend .depend: $(SRCS) $(MODS)
	makedepf90 -o abc_db.exe  $(SRCS) $(MODS)  > .depend

# To make "forced_compile" the default target, ie work with just "make",
# move it direcely below the "OBJS" definition and
# comment out the "include .depend" and "depend" targets
.PHONY: forced_compile
forced_compile:
	@(echo "Forced compile of all routines")
	$(FC) -c $(FFLAGS) $(MODS)
	$(FC) -c $(FFLAGS) $(SRCS)
	$(FC) $(OBJS) -o $(EXE)

.PHONY: clean
clean:
	@(echo "Cleaning abc directory ...")
	@(echo "removing *.o *.mod and *.exe files")
	-rm -f *.o
	-rm -f *.mod
	-rm -f *.exe
	-rm -f .depend
	

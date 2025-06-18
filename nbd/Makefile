##############################################################
###                                                        ###
###           Make File for the NBD Code                   ###
###                                                        ###
##############################################################

.KEEP_STATE:

.SUFFIXES: .f90 .o
FC = gfortran 
FFLAGS = -O3 -ffree-line-length-none -C

FOBS = global_parameters.o global_variables.o timing_info.o accel_loop_list.o \
	nbd.o allocate_memory.o \
	get_accel_info_faster.o \
	get_accel_for_active_faster.o \
	get_ics.o \
	hermite_4thorder_PECn_fixed_dt.o \
        hermite_4thorder_PECn_block_dt.o \
	hermite_4thorder_PECn_global_dt.o \
        verlet_2ndorder_global_dt.o \
	initialise_run.o \
	read_param_file.o write_a_snapshot_dir.o \
	diagnostics.o get_energies_and_momentum.o \
	synchronise_posvel_for_output.o

EXEC = ~/home/dxs_cosmos/bin/nbd

.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

$(EXEC): $(FOBS)
	$(FC) $(FOBS) $(FFLAGS) -o $(EXEC)

clean:
	\rm -f *.o *.mod $(EXEC)

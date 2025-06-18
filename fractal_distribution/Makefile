.KEEP_STATE:

.SUFFIXES: .f90 .o
FC = gfortran 
FFLAGS = -O3 -ffixed-line-length-132 -C

FOBS = global_parameters.o global_variables.o \
	whole_code.o  \
	ICS_test.o

EXEC = ~/bin5/gen_ics

.f90.o :
	$(FC) $(FLAGS) -c $*.f90

$(EXEC): $(FOBS)
	$(FC) $(FOBS) $(FFLAGS) -o $(EXEC)

clean:
	\rm -f *.o *.mod $(EXEC)
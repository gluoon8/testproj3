
# Fortran compiler
FC=gfortran

# Code location
VPATH=src/

######################
## 	 HELP	    ##
######################
##
.PHONY. : help
help:
	@grep "^##" Makefile

# ~ LINKING ~
## 	· make: Compiles the program.
##

all: MD.x
MD.x:  pbc.o integrators.o forces.o initialize.o  main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@


# ~ RUNNING THE PROGRAM ~
## 	· make run: Runs the program. If the program is not compiled, 
##			  it compiles it first.
##

.PHONY.: run
run: MD.x
	./MD.x

# ~ COMPILING ~

pbc.o: pbc/pbc.f90
	$(FC) -c $^

integrators.o: integrators/integrators.f90
	$(FC) -c $^

forces.o: forces/forces.f90
	$(FC) -c $^

initialize.o: initialization/initialize.f90
	$(FC) -c $^

main.o: main.f90
	$(FC) -c $^



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
## 	· make comp: Compiles the program.
##
.PHONY.: comp
comp: MD.x
MD.x:  pbc.o integrators.o forces.o initialize.o  main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@


# ~ RUNNING THE PROGRAM ~
## 	· make run: Runs the program. If the program is not compiled, 
##			  it compiles it first.
##

.PHONY.: run
run: MD.x
	./MD.x


# ~ CLEANING OUTPUT FILES ~
## 	· make out: Removes ONLY the output files.
##
.PHONY.: out
out:
	rm -f *.dat

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

# ~ PLOTS ~
## 	· make plot: Plots the output data:
##		- Epot, Ekin, Etot vs time
##		- Momentum vs time
##		- T , pressure vs time

.PHONY.: plot
plot:
	@python3 ./src/plots/plots.py

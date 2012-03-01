# PREDEFINED VARIABLES
FFLAGS = -c -Wall
FORT = mpif90

# MODULES 
MODULES  = param arrays
MODULESF90  = $(MODULES:=.f90)
MODULESO    = $(MODULES:=.o)

#ROUTINES
ROUTINES = parallel_initial_data allocate initialize parallel_grid grid initialdata observables output_0D parallel_output_1D output_1D store_levels sources internal_boundaries external_boundaries evolution step_rk4 
ROUTINESF90 = $(ROUTINES:=.f90)
ROUTINESO   = $(ROUTINES:=.o)

#MAIN PROGRAM
PROGRAM  = main

#EXECUTABLE
EXEC = exec.out


# LINKING
$(PROGRAM).out: $(MODULES) $(ROUTINES) $(PROGRAM).o
	$(FORT) $(PROGRAM).o $(MODULESO) $(ROUTINESO) -o $(EXEC)
	rm -rf $(MODULESO) $(ROUTINESO) $(PROGRAM).o $(MODULES).mod

# MODULES' RULES
$(MODULES): $(MODULESF90) 
	$(FORT) $(FFLAGS) $@.f90 -o $@.o

# ROUTINES' RULES
$(ROUTINES): $(ROUTINESF90) 
	$(FORT) $(FFLAGS) $@.f90 -o $@.o

# PROGRAM COMPILATION
$(PROGRAM).o: $(PROGRAM).f90 
	$(FORT) $(FFLAGS) $< -o $@

clean:
	rm -rf *.out *.o *.f90# *.log *.mod *.f90~

run:
	mpirun -np 4 ./$(EXEC)

all: clean $(PROGRAM).out run
	gnuplot script


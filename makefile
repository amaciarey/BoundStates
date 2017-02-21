#f90=ifort -O3
f90=gfortran -Wall -O3

bound_states:bound_states.o Shoot.o Bisection.o
	$(f90) -o bound_states *.o
bound_states.o:bound_states.f90
	$(f90) -c bound_states.f90
Shoot.o:Shoot.f90
	$(f90) -c Shoot.f90
Bisection.o:Bisection.f90
	$(f90) -c Bisection.f90

clean:
	rm -f *.o

program bound_states

implicit none

real (kind=8)    :: Potential
real (kind=8)    :: xmax,tol,d
integer (kind=4) :: N,Ne,itmax,m

real (kind=8)    :: E0,E1,E,Emin,Emax,Etrial,dE
real (kind=8)    :: Delta_0,Delta_1
real (kind=8)    :: x,h
integer (kind=4) :: i
logical          :: bound

real (kind=8),dimension (:),allocatable :: V,u

!Reading input parameters

open (unit=1,file='bound_states.in',status='old')

read (1,*) xmax
read (1,*) N
read (1,*) tol
read (1,*) itmax
read (1,*) m
read (1,*) d
read (1,*) Emin
read (1,*) Emax
read (1,*) Ne

close (unit=1)

!Definition of aditional parameters

h  = xmax/real(N)

!Definition of the array that contain the tabulated potential

allocate (V(N),u(N))

open (unit=4,file='potential.out')

do i=1,N
   x    = real(i)*h
   V(i) = Potential(m,d,x)
   write (4,*) x,V(i)
end do

close (unit=4)

!This Emin is only valid for the bilayer dipole-dipole interaction

dE = (Emax-Emin)/real(Ne-1)

open (unit=2,file='eigenfunctions.out')
open (unit=3,file='eigenvalues.out')

Etrial = Emin

do while (Etrial<Emax)

   E0 = Etrial
   E1 = Etrial+dE

   !Integration of the differential equation with E = E0

   call Shoot(N,h,E0,V,u,Delta_0)

   !Integration of the differential equation with E = E1

   call Shoot(N,h,E1,V,u,Delta_1)

   !Bisection method to find an eigenvalue of the energy in the 
   !range [E0,E1]

   call Bisection(N,h,tol,itmax,Delta_0,Delta_1,E0,E1,V,u,E,bound)

   !Outputs

   if (bound) then
      
      write (3,*) 'E=',E
     
      do i=1,N
         x = real(i)*h
         write (2,'(2g20.10e3)') x,u(i)
      end do
      write (2,*)
      write (2,*)

   end if

   Etrial = Etrial+dE
  
end do

close (unit=2)
close (unit=3)

deallocate (V,u)

end program bound_states

!=======================================================================
!      ROUTINES
!=======================================================================

real (kind=8) function Potential(m,h,x)

implicit none

real (kind=8)    :: x,h
integer (kind=4) :: m

!Hydrogen atom

Potential = -1.d0/abs(x)+real(m*(m+1))/(x*x)

!Potential = real(m*(m+1))/(x*x)+4.d0*((1.d0/x)**12-(1.d0/x)**6)

!Bilayer dipole-dipole system

!Potential = (real(m*m)-0.25d0)/(x*x)+(x*x-2.d0*h*h)/sqrt(x*x+h*h)**5

!Two dimensional square well

!Potential = (real(m*m)-0.25d0)/(x*x)

!if (x<=1.d0) Potential = Potential-h

!3D square well potential

!Potential = real(m*(m+1))/(x*x)

!if (x<=1.d0) Potential = Potential-5.d0

return
end function Potential

!-----------------------------------------------------------------------

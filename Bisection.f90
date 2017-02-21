!======================================================================
!         SUBROUTINE: BISECTION
!======================================================================
!
!This subroutine solves the equation:
!
!         BETA_right(E) = BETA_left(E)
!
!using the bisection method. Here BETA is the logarithmic derivative 
!of the wave function at the matching point.
!The subroutine receives two values of the energy, E0 and E1, and
!checks if there are a root of the equation between E0 and E1, if 
!there are no roots, the routine does nothing. If there are any root
!it is calculated using the bisection method up to a given precision.
!The subroutine also normalizes the wave function given by Shoot 
!using the trapezoidal rule. 
!
!The arguments of the subroutine are the following:
!
!   > N       : Number of points in the integration grid
!   > h       : Stepsize of the integration
!   > tol     : Tolerance required in the precision of the eigenenergy
!   > itmax   : Maximum number of bisections to achieve tol
!   > Delta_0 : BETA_right(E0)-BETA_left(E0)
!   > Delta_1 : BETA_right(E1)-BETA_left(E1)
!   > E0,E1   : Initial values of energy to perform the bisection
!   > V       : Array that contains the values of the effective 
!               potential evaluated on the grid points. Dim(V) = N
!   > u       : Array that contains the values of the radial wave
!               function evaluated on the grid points. Dim(u) = N
!   > E       : Eigenvalue of the Schrodinger equation calculated 
!               using the bisection method
!   > bound   : Logical variable that tell the main program if 
!               a bound state was found, bound = .true. implies
!               that there are a bound state while bound = .false.
!               implies that there are no bound states between 
!               E0 and E1.
!
!-----------------------------------------------------------------------


subroutine Bisection(N,h,tol,itmax,Delta_0,Delta_1,E0,E1,V,u,E,bound)

implicit none

real (kind=8)    :: h,tol,Delta_0,Delta_1,E0,E1,E
integer (kind=4) :: N,itmax

real (kind=8)    :: E2,Delta_2,Error,norm
integer (kind=4) :: istep,i
logical          :: bound

real (kind=8),dimension (N) :: V,u

!Check if there are any bound state between E0 and E1

if (Delta_0*Delta_1>0.d0) then

   !If there are no bound state exit the subroutine

   bound = .False.

else if (Delta_0*Delta_1<0.d0) then

   !If there are a bound state divide the domain [E0,E1] in two 
   !parts defining E2 = (E0+E1)/2 and look for a bound state in 
   !the two new subdomains: [E0,E2],[E2,E1]

   Error = abs(Delta_0-Delta_1)
   istep = 1

   bound = .False.

   do 
      
      E2 = 0.5d0*(E0+E1)
      
      !Define the new function q for the energy E2

      call Shoot(N,h,E2,V,u,Delta_2)

      if (Delta_2*Delta_0<0.d0) then
         E1      = E2
         Delta_1 = Delta_2
         Error   = abs(Delta_1-Delta_0)
      else
         E0      = E2
         Delta_0 = Delta_2
         Error   = abs(Delta_1-Delta_0)
      end if

      istep = istep+1

      if (Error<tol) exit
      if (istep>=itmax) exit
     
   end do 

   if (Error<tol) then
      bound = .True.
      E     = E2
   else
      bound = .False.
   end if

   !Normalizing the obtained wave function

   norm = 0.d0

   do i=1,N
      norm = norm+u(i)*u(i)
   end do

   norm = sqrt(h*norm)

   do i=1,N
      u(i) = u(i)/norm
   end do

end if

!Print results of the bisection calculation

if (bound) then
   print *, '-----------------------------------------------'
   print *, '# E          =',E
   print *, '# Bisections =',istep
   print *, '# Delta      =',Error
end if

return 
end subroutine Bisection

!-----------------------------------------------------------------------

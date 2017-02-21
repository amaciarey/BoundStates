!=======================================================================
!         SUBROUTINE: SHOOT
!=======================================================================
!
!This subroutine solves the Schrodinger equation using the Numerov's 
!algorithm and the shooting method. The Schrodinger equation for a bound
!state is an eigenvalue equation, this means that not only the wave 
!function u(x) must be evaluated, also the eigenvalue E. Then the 
!solution of the problem is given by the pair {E,u(x)}. This subroutine
!takes a trial value for the eigenenergy and then integrates the radial
!Schrodinger equation for this energy. The integration is carried by
!splitting the space in two different regions, (i) From xmin to xmatch, 
!and (ii) From xmatch to xmax, where xmatch is the classical right 
!turning point.
!The subroutine uses the Numerov's method to perform the integration in 
!the two regions:
!
!   (i)  From xmin to xmatch integrates forward with a step h.
!   (ii) From xmax to xmatch integrates backward with a step -h
!
!Then continuity is imposed at xmatch and finally the logarithmic 
!derivatives of the in and out solutions is evaluated.
!
!The arguments of the subroutine are:
!
!   > N     : Number of points in the integration grid
!   > h     : Integration step
!   > E     : Trial value for the energy
!   > V     : Array that contains the values of the effective potential
!             in the grid points. dim(V) = N
!   > u     : Array that contains the values of the radial wave function
!             in the grid points. dim(u) = N
!   > Delta : Difference between the logarithmic derivative of the wave
!             function from the forward and backward integration.
!
!-----------------------------------------------------------------------

subroutine Shoot(N,h,E,V,u,Delta)

implicit none

real (kind=8)    :: h,E,Delta
integer (kind=4) :: N,im,i

real (kind=8)    :: step,eps,gafter,gpresent,gbefore
real (kind=8)    :: left_der,right_der,scale,um
real (kind=8)    :: a11,a12,a21,a22,b1,b2,dudx

real (kind=8),dimension (N) :: q,u,V

step  = h*h/12.d0
eps   = 1.d-3
im    = 1

!Define the function q(x) = V(x)-E and look for the classical
!right turning point

do i=1,N
   q(i) = (V(i)-E)
end do

do i=1,N-1
   if ((q(i)*q(i+1)<0.d0).and.(V(i)<E)) im = i
end do

!Define the conditions on the boundaries of the integration grid

u(1) = h
u(2) = 2.d0*h

u(N)   = 0.d0
u(N-1) = eps

!Integrate from left to right until im is reached

do i=3,im
   gafter   = 1.d0-step*q(i)
   gpresent = 2.d0*(1.d0+5.d0*step*q(i-1))
   gbefore  = 1.d0-step*q(i-2)
   u(i)     = gpresent*u(i-1)-gbefore*u(i-2)
   u(i)     = u(i)/gafter
end do

um = u(im)

!Evaluation of the derivative at the matching point from the left

left_der = (u(im)-u(im-1))/h

!Integrate from right to left until im is reached

do i=N-2,im,-1
   gafter   = 1.d0-step*q(i)
   gpresent = 2.d0*(1.d0+5.d0*step*q(i+1))
   gbefore  = 1.d0-step*q(i+2)
   u(i)     = gpresent*u(i+1)-gbefore*u(i+2)
   u(i)     = u(i)/gafter
end do

!Scale the right solution to have a continuous wave function

scale = um/u(im)

do i=im,N
   u(i) = scale*u(i)
end do

!Evaluation of the derivative at the matching point from the right 

right_der = (u(im+1)-u(im))/h

!Difference of the right and left derivatives 

Delta = right_der-left_der

return
end subroutine Shoot

!-----------------------------------------------------------------------

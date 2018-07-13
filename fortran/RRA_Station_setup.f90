! EXAMPLE OF SET UP TEST CASE 
! by Letizia Tedesco (letizia.tedesco@ymparisto.fi) 
! January 2007, Revised November 2012 
! Variable declarations 
real(MK) :: dt  !<  
real(MK) :: Fwater  !<  
real(MK) :: hi_min  !<  
real(MK) :: hmi_min  !<  
real(MK) :: hs_min  !<  
real(MK) :: hs_prec_min  !<  
integer :: nmax  !<  
real(MK) :: Tfreez  !<  
real(MK) :: ZERO  !<  
! 

!	EXAMPLE OF AN ARTIC SIMULATION: 27/11/2005 00:00 to 27/06/2006 65.10 N, 308.32 E 
! Forcing data interpolated by Gradsnc on the location of station 
!	and interpolated by MATlab with a linear interpolation 

!!!!!!!!!!!!!!!!!!!!!!!!!!! !!! TO BE SET UP BY THE USER !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
dt=5400 ! = deltat = time step (s)
nmax=5837 ! maximum number of time steps (-)
hs_prec_min=0.02 ! minimum snow precipitation (m) to activate the model !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hs_min=0.02 ! minimum snow thickness (m) to activate the model !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hmi_min=0.05 ! minimum snow ice/supeimposed ice thickness (m) to activate the model !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hi_min=0.05 ! minimum sea ice thickness (m) to activate the model !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
ZERO=10**-5 ! ZERO OF THE MODEL !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
Fwater=9.3 ! = Fw = oceanic heat fluxes (W m-2) !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS PARAMETER!!!
Tfreez=271.3722 ! = Tfr = -0.4 C seawater freezing temperature (K) !!!BE AWARE: THIS STRONGLY DEPENDS ON THE LOCATION!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

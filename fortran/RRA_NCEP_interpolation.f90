! EAXMPLE OF INPUT DATA INTERPOLATION FROM 6H TO 1.5 H 
! by Letizia Tedesco (letizia.tedesco@ymparisto.fi) 
! January 2007, Revised November 2012 
! !!! BE AWARE: TO BE RUN BEFORE MODEL RUN !!! 
! Variable declarations 
real(MK), dimension(1) :: Cl  !< !m2f: check dim(:) 
real(MK), dimension(1) :: Fsd_cloud  !< !m2f: check dim(:) 
real(MK), dimension(1) :: P_rate  !< !m2f: check dim(:) 
real(MK), dimension(1) :: qa  !< !m2f: check dim(:) 
real(MK), dimension(1) :: qs  !< !m2f: check dim(:) 
real(MK), dimension(1) :: Ta  !< !m2f: check dim(:) 
real(MK), dimension(1) :: U  !< !m2f: check dim(:) 
real(MK) :: Ua  !<  
real(MK), dimension(1) :: V  !< !m2f: check dim(:) 
real(MK) :: x  !<  
real(MK) :: xi  !<  
! 

! !!! BE AWARE: PATH TO FORCING TO BE SET BY THE USER !!! !!!!!!!!!!!!!!!!!!!!!!!!!!! 
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/Cl.txt !(total cloud cover in percentage)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/Fsd_cloud.txt !(downwars shortwave)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/P_rate.txt !(precipitation rate)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/qa.txt !(specific humidity at 2 m height)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/qs.txt !(specific humidity at the surface)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/Ta.txt !(air temperature)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/U.txt !(U wind speed at 10 m height)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/3/V.txt !(V wind speed at 10 m height)

x=1:4:5840 !!! to be set up by user
xi=1:1:5840 !!! to be set up by user

!m2f: Cl=Cl'
Cl=transpose(Cl)
Cl=interp1(x,Cl,xi,'linear')
!m2f: Cl=Cl'/100 !!! convert cloudiness from percentage to fraction
Cl=transpose(Cl)/100 !!! convert cloudiness from percentage to fraction
Cl(5838:5840)=[ ] 

!m2f: Fsd_cloud=Fsd_cloud'
Fsd_cloud=transpose(Fsd_cloud)
Fsd_cloud=interp1(x,Fsd_cloud,xi,'linear')
!m2f: Fsd_cloud=Fsd_cloud'
Fsd_cloud=transpose(Fsd_cloud)
Fsd_cloud(5838:5840)=[ ] 

!m2f: qa=qa'
qa=transpose(qa)
qa=interp1(x,qa,xi,'linear')
!m2f: qa=qa'
qa=transpose(qa)
qa(5838:5840)=[ ] 

!m2f: qs=qs'
qs=transpose(qs)
qs=interp1(x,qs,xi,'linear')
!m2f: qs=qs'
qs=transpose(qs)
qs(5838:5840)=[ ] 

!m2f: Ta=Ta'
Ta=transpose(Ta)
Ta=interp1(x,Ta,xi,'linear')
!m2f: Ta=Ta'
Ta=transpose(Ta)
Ta(5838:5840)=[ ] 

!m2f: P_rate=P_rate'
P_rate=transpose(P_rate)
P_rate=interp1(x,P_rate,xi,'linear')
!m2f: P_rate=P_rate'
P_rate=transpose(P_rate)
P_rate(5838:5840)=[ ] 

!m2f: U=U'
U=transpose(U)
U=interp1(x,U,xi,'linear')
!m2f: U=U'
U=transpose(U)
U(5838:5840)=[ ] 
!m2f: V=V'
V=transpose(V)
V=interp1(x,V,xi,'linear')
!m2f: V=V'
V=transpose(V)
V(5838:5840)=[ ] 
Ua=sqrt(U.**2+V.**2)

clear U V x xi

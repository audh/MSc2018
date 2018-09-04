%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERA INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   !!! PATH TO ERA FORCING FILES TO BE SET BY USER !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/Cl.txt           %(total cloud cover in percentage)
load /home/riesna/MSc2018/ESIM_AntTest/NCEP/1/Fsd_cloud.txt    %(downwars shortwave)
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/P_rate.txt       %(precipitation rate)
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/qa.txt           %(specific humidity at 2 m height)
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/qs.txt           %(specific humidity at the surface)
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/Ta.txt           %(air temperature)
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/U.txt            %(U wind speed at 10 m height)
load /home/riesna/MSc2018/ESIM_AntTest/ERA/1/V.txt            %(V wind speed at 10 m height)

x=1:4:5840; %%% to be set up by user
xi=1:1:5840; %%% to be set up by user

Cl=Cl';
Cl=interp1(x,Cl,xi,'linear');
Cl=Cl'/100; %%% convert cloudiness from percentage to fraction
Cl(5838:5840)=[]; 

Fsd_cloud=Fsd_cloud';
Fsd_cloud=interp1(x,Fsd_cloud,xi,'linear');
Fsd_cloud=Fsd_cloud';
Fsd_cloud(5838:5840)=[];

qa=qa';
qa=interp1(x,qa,xi,'linear');
qa=qa';
qa(5838:5840)=[];

qs=qs';
qs=interp1(x,qs,xi,'linear');
qs=qs';
qs(5838:5840)=[];

Ta=Ta';
Ta=interp1(x,Ta,xi,'linear');
Ta=Ta';
Ta(5838:5840)=[];

P_rate=P_rate';
P_rate=interp1(x,P_rate,xi,'linear');
P_rate=P_rate';
P_rate(5838:5840)=[];

U=U';
U=interp1(x,U,xi,'linear');
U=U';
U(5838:5840)=[];
V=V';
V=interp1(x,V,xi,'linear');
V=V';
V(5838:5840)=[];
Ua=sqrt(U.^2+V.^2);

clear U V x xi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
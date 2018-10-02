%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESIM3: STANDARD RUN FOR ANTARCTICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riesna R Audh

% Revised script from Letizia Tedesco (2007)

% Standard run of the ESIM2 for the Antarctic to produce the BFM-SI input file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATION SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  !!! TO BE SET UP BY THE USER !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=5400;                               % = deltat = time step (s)
nmax=5837;                             % maximum number of time steps (-)
hs_prec_min=0.02;                      % minimum snow precipitation (m) to activate the model                   !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hs_min=0.02;                           % minimum snow thickness (m) to activate the model                       !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hmi_min=0.05;                          % minimum snow ice/superimposed ice thickness (m) to activate the model   !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hi_min=0.05;                           % minimum sea ice thickness (m)  to activate the model                   !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
ZERO=10^-5;                            % ZERO OF THE MODEL                                                      !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
Fwater=8.25;                            % = Fw = oceanic heat fluxes (W m-2)                                     !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS PARAMETER!!!
Tfreez=271.35;                         % = Tfr = -1.8 C seawater freezing temperature (K)                       !!!BE AWARE: THIS STRONGLY DEPENDS ON THE LOCATION!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL PARAMETERS AND CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL PROGNOSTIC VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   deltahi_bott        (m)             variation in thickness due to sea ice melting at the bottom
%   deltahi_melt_surf   (m)             variation in thickness due to sea ice melting at the surface
%   deltahmi_melt_surf  (m)             variation in thickness due to snow ice/supeimposed ice melting at the surface
%   deltahs_melt_surf   (m)             variation in thickness due to snow melting at the surface
%	hs                  (m)             1 layer of new snow + 1 layer old snow
%	hmi                 (m)             1 layer of snow ice + 1 layer of superimposed ice
%	hi                  (m)             1 layer of sea ice
%   T0                  (K)             surface temperature
%   Sbr                 (per mill)      brines salinity
%	Sice                (per mill)      sea ice bulk salinity
%   Vbr_ice             (non-dim)       sea ice brines volume fraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL FORCING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Cl          (non-dim)   total cloud cover fraction
%	Fsd_cloud   (W m-2)     downward solar radiation flux
%	P_rate      (m s-1)     precipitation rate
%   qa          (non-dim)   specific humidity of the air
%   qs          (non-dim)   specific humidity of the surface
%	Ta          (K)         air temperature at 10 m height
%   Ua          (m s-1)     wind speed at 10 m height

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL DIAGNOSTIC VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	alpha           (non-dim)       albedo
%   cpa             (J kg^-1 K^-1)  sea ice heat capacity
%   delta_bucket    (m)             variations in snow in the bucket model
%   ea              (non-dim)       air relative humidity
%   Fbott           (W m-2)         conductive flux at the bottom
%   Fl              (W m-2)         net longwave radiation flux
%   Fla             (W m-2)         net latent heat flux
%   Fld             (W m-2)         downward longwave radiation flux
%   Flu             (W m-2)         upward longwave radiation flux
%   Fs              (W m-2)         net solar radiation flux
%   Fse             (W m-2)         net sensible heat flux
%   Fsurf           (W m-2)         conductive flux at the surface
%   Fw              (W m-2)         oceanic heat flux
%	hs_prec         (m)             snow precipitation
%   hs_prec_bucket  (m)             snow precipitation in the bucket model
%   ice_fr          (non-dim)       ice fraction: 0.0 if open ocean, 1.0 if sea ice cover
%   I0              (W m-2)         fraction of solar radiation penetrating the surface
%   IM              (W m-2)         fraction of solar radiation penetrating snow ice/supeimposed ice
%   ISI             (W m-2)         fraction of solar radiation penetrating sea ice below the first 10 cm
%   ISI_10          (W m-2)         fraction of solar radiation penetrating sea ice in the above first 10 cm
%   ISI_layer       (W m-2)         fraction of solar radiation penetrating the biologically-active layer sea ice in the above first 10 cm
%   ISI_bio         (W m-2)         fraction of solar radiation penetrating sea ice and reaching the middle point of the biologically active layer
%   ki              (W m^-1 s^-1)   thermal conductivity of sea ice
%   kmi_clear       (m-1)           snow ice/superimposed ice extiction coeff in clear sky conditions
%   kmi_cloudy      (m-1)           snow ice/superimposed ice extiction coeff in cloudy sky conditions
%   kmi_av          (m-1)           mean snow ice/superimposed ice extiction coeff as function of Cl
%   ks              (W m^-1 s^-1)   thermal conductivity of snow
%   ks_clear        (m-1)           snow extiction coeff in clear sky conditions
%   ks_cloudy       (m-1)           snow extiction coeff in cloudy sky conditions
%   ks_snow           (m-1)           mean snow extiction coeff as function of Cl
%   ksi_10_clear    (m-1)           sea ice extiction coeff at 10 cm in clear sky conditions
%   ksi_10_cloudy   (m-1)           sea ice extiction coeff at 10 cm in cloudy sky conditions
%   ksi_10_av       (m-1)           mean sea ice extiction coeff at 10 cm as function of Cl
%   ksi_clear       (m-1)           sea ice extiction coeff in clear sky conditions
%   ksi_cloudy      (m-1)           sea ice extiction coeff in cloudy sky conditions
%   ksi_av          (m-1)           mean sea ice extiction coeff as function of Cl
%   roi_surf        (kg m-3)        density of pure ice at the surface
%   ro_sice_surf    (kg m-3)        bulk density of sea ice at the surface as function of T,S;
%   ro_sice_bott    (kg m-3)        bulk density of sea ice at the bottm (density of pure ice at the bottom is fixed at 0.900)
%   roi_av          (kg m-3         density of pure ice of the layer as function of Tice_av in g/cm^3
%   ro_sice_bulk    (kg m-3)        bulk density of ice as function of Tice,Sice in g/cm^3
%   ro_br           (kg m-3)        brines density as function of brines salinity in g/cm^3
%   ro_br_5         (kg m-3)        brines density in the biologically active layer as function of brines salinity in g/cm^3
%   ro_br_bio       (kg m-3)        mean brines density in the biologically active layer as function of brines salinity in g/cm^3
%   T_beg           (K)             initial guess for temperature iteration
%   T_iter          (K)             iteration temperature
%   Tice            (??C)            temperature of sea ice at the upper interface
%   Tice_5          (??C)            temperature of sea ice at brines volume 0.05
%   Tice_bio        (??C)            temperature of sea ice in the biological active layer
%   Tis             (K)             temperature at the interf snow ice (or supeimposed ice)/sea ice
%   Tss             (K)             temperature at the interf snow/snow ice (or supeimposed ice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	SET PARAMETRES AND CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_ow=0.15;          %   (non-dim)       seawater albedo
alpha_mi=0.55;          %   (non-dim)       intermediate layer albedo
beta=2.0;%2.667;		%   (non-dim)       empir coeff from snow ice to ice (after Lepparanta,1983)
beta2=0.13;             %   (W m-1)         emp cost (after Untersteiner, 1961)
c0=2093.0; 				%   (J (kg*??deg)-1) specific heat of fresh ice
c1=21.875;				%   (K)             emp.costant
c2=7.66;				%   (K)             emp.costant
c10=0.1;                %   (m)             (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
c11=0.44;               %   (m^-0.28)       (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
c12=0.075;              %   (m^-2)          (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
cail=1.75*10^-3;		%   (non-dim)       bulk transfer coefficient
cais=1.20*10^-3;		%   (non-dim)       bulk transfer coefficient
cpair=1004.0;           %   (J kg^-1 K^-1)  air heat capacity
cpw=4186.0;             %   (J kg^-1 K^-1)  specific heat of sea water
cs=2093.0;				%   (J (kg K)-1)    specific heat of snow
cwih=2.8*10^-4;			%   (non-dim)       bulk heat trans coeff (after Omstedt and Wettlaufer,1992)
deltat=dt;              %   (s)             time step 1.5 h
emis=0.97;				%   (non-dim)       water/ice emissivity
emp=0.62197;            %   (mbar)          empirical coefficient for specific humidity computation
gamm=18.0;              %   (J (K*kg)-1)    emp coeff (after Bitz and Lipscomb, 1999; otherwise 17.2 after Untersteiner, 1961)
infra_fr=0.77;          %   (non-dim)       Infrared partition as for water type II (after Jerlov, 1968)
k0=2.03;                %   (W m^-1 K^-1)   fresh ice thermal conductivity
ks_prec=0.0641;         %   (W m^-1 K^-1)   precipit thermal conductivity 
kw=0.563;               %   (W m^-1 K^-1)   sea water thermal conductivity (after Lepparanta, 1983)
kmi=0.9;				%   (W m^-1 K^-1)   snow ice thermal conductivity (after Lepparanta, 1983)
k_ocean_red=0.6667;     %   (m^-1)          seawater extiction coefficient (infrared, water type II, after Jerlov 1968)
k_ocean_vis=0.0714;     %   (m^-1)          seawater extiction coefficient (visible, water type II, after Jerlov 1968)
h_mix=10.0;             %   (m)             mixed layer depth
Lai=2.501*10^6;			%   (J kg^-1)       latent heat of vaporization
L0i=297000.0;           %   ()              latent heat of fusion of pure ice
L0m=315000.0;           %   ()              latent heat of fusion of meteoric ice
L0s=334000.0;           %   ()              latent heat of fusion of pure snow
Lv=2.501*10^6;          %   (J kg-1)        latent heat of vaporization of fresh water at 273.15 K
mu=0.054; 				%   (??C)            ratio between Tfr and brines salinity (after Assur, 1958)
P0=1013.0;				%   (mbar)          seawater pressure at the surface
q1=1.16378*10^7;        %   (kg m-3)        emp coeff
q2=5897.8;              %   (K)            emp coeff
roa=1.225;              %   (kg m-3)        air density
roi_bott=0.9;           %   (Kg m-3)        pure sea ice density at the bottom
romi=850.0;             %   (kg m-3)        snow ice/superimposed ice density
rosa=1.68*10^-5;        %   1.68*10^-7;    %5.88*10^-8;        %   (C??^-1*m*s^-1)  emp coeff (after Vancoppenolee et al, 2007; otherwise 1.68*10^-7 after Cox and Weeks, 1988)
ros_prec=250;           %   (kg m-3)        snow precip density
row=1026.0;             %   (kg m-3)        water density
Sbr_end=18.5185;         %  (Kg m-3)        pure brines density at the bottom
sigma=5.68*10^-8;		%   (W m^-2 K)      Stefan-Boltzmann constant
si_fract=0.0;           %   (non-dim)       fraction of snow transformation in snow ice !!!IT CAN BE SET UP BY THE USER DEPENDING ON SNOW TYPE, ECC.)!!!
ss_fract=0.0;           %   (non-dim)       fraction of snow transformation in superimposed ice (after Cheng et al, 2006)!!!IT CAN BE SET UP BY THE USER DEPENDING ON SNOW TYPE, ECC.)!!!
%Sw=32.0;				%   (per mill)      water salinity
Ssnow=0.0;              %   (per mill)      snow salinity
Tb=273.155;				%   (K)             empirical temp
Tfr=Tfreez;             %   (K)             freezing temp of sea ice
Tfrs=273.15;			%   (K)             freezing temp of snow/fresh water
Va=0.015;               %   ()              gas volume fraction in the ice is fixed
viola=20.0;             %                   for salinity compuatation (after  Vancoppenolee et al, 2007)
vis_fr=0.7;             %                   fraction of visible light penetrating the surface layer when SIM is on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!!TO BE SET BY THE USER!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmix(1)=274.55;                     % Initial mixed layer temperature
T0(1)=Tmix(1);                      % Initial surface temperature
Ts(1)=Tfr;                          % Initial snow temperature
Tsi(1)=Tfr;                         % Initial snow ice temperature
Ti(1)=Tfr;                          % Initial sea ice temperature

hi(1)=ZERO;                         % Initial sea ice thickness
hs(1)=ZERO;                         % Initial snow thickness
hmi(1)=ZERO;                        % Initial snow ice/superimposed ice thickness
hmi_new(1)=0.0;
hs_prec_bucket(1)=0.0;              % Initial snow thickness in the bucket model
hs_prec(1)=0.0;                     % Initial snow precipitation
hs_tot(1)=0.0;
delta_sublim(1)=0.0;
delta_bucket(1)=0.0;

snow_fr(1)=0.0;
sea_water_fr(1)=1.0;

ki(1)=1.90;                         % Initial sea ice thermal conductivity
ks(1)=0.25;                         % Initial snow thermal conductivity
ksi_10_av(1)=-6.6*Cl(1)+17.1;       % Initial mean of sea ice extiction coefficient (10 cm layer)
ks_av(1)=0.20;
cpi(1)=2200.0;                      % Initail sea ice heat capacity
ros(1)=250;                         % Initial snow density
ros_new(1)=150;                     % Initial new snow density

Vbr_ice(1)=0.0;                     % Initial brines volume
Sbr_ice(1)=0.0;                     % Initial brines salinity
Sice(1)=0.0;                        % Initial sea ice bulk salinity
Tice(1)=Tfr-273.15;                 % Initial sea ice temperature
Sice_bott(1)=0.0;
Ssnowice(1)=0.0;
Vbr_i(1)=0.0;                       % Initial internal brines volume
Sbr_i(1)=0.0;                       % Initial internal brines salinity
Si(1)=0.0;                          % Initial internal sea ice bulk salinit

Sice_5(1)=0.0;                      % Initial sea ice bulk salinity at Vbr=0.05
Sbr_5(1)=0.0;                       % Initial brines salinity at Vbr=0.05
Tice_5(1)=Tfr-273.15;                   % Initial sea ice temperature at Vbr=0.05
hi_5(1)=ZERO;                       % Initial sea ice thickness at Vbr=0.05

Sice_bio(1)=0.0;                    % Initial sea ice bulk salinity in the bio layer
Tice_bio(1)=0.0;                    % Initial sea ice temperature in the bio layer
Sbr_bio(1)=0.0;                     % Initial brines salinity in the bio layer
hi_bio(1)=hi_5(1);                  % Initial sea ice thickness in the bio layer
Vbr_bio(1)=0.0;                     % Initial brines volume in the bio layer
ISI_bio(1)=0.0;                     % Fraction of solar radiation penetrating sea ice and reaching the middle point of the biologically active layer

R(1)=0.0;

ro_sice_bott(1)=0.917;              % Initial sea ice bulk density at the bottom
ro_sice_surf(1)=0.900;              % Initial sea ice bulk density at the surface
ro_sice_bulk(1)=0.910;

ks_snow(1)=25.0;                    % Initial mean of snow extiction coefficient
%kmi_av(1)=-3.3*Cl(1)+21.05;         % Initial mean of snow ice/superimposed ice extiction coefficient
%ksi_10_av(1)=-6.6*Cl(1)+17.1;       % Initial mean of sea ice extiction coefficient (10 cm layer)
%ksi_av(1)=-0.1*Cl(1)+1.6;         % Initial mean of sea ice extiction coefficient
alpha(1)=0.15;                      % Initial albedo value

ki_5_bott(1)=2.0;                   % Initial sea ice thermal conductivity between Vbr=0.05 and the bottom
ki_ice_5(1)=2.0;                    % Initial sea ice thermal conductivity between sea ice surface and Vbr=0.05
ki_ice_bio(1)=2.0;
ki_bio_bott(1)=2.0;                 % Initial sea ice thermal conductivity between the bio layer and the bottom
ki_5_bio(1)=2.0;                    % Initial sea ice thermal conductivity between Vbr=0.05 and the bio layer
ki_ice_bott(1)=2.0;                 % Initial sea ice thermal conductivity at the bottom

Sw(1)=33.6;                         % Seawater salinity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M O D E L S T A R T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=2:nmax                % maximum number of timestep
    Sw(i)=33.6;                     % Seawater salinity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   SNOW DENSITY
    
    if T0(i-1)<Tfr
        ros(i)=300;
        ros_new(i)=250;
    else
        ros(i)=350;
        ros_new(i)=300;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %  SNOW PRECIPITATION
    if (Ta(i)<Tfr)
        if hi(i-1)>hi_min
            hs_prec(i)=(P_rate(i)/ros_prec)*deltat;
        else
            hs_prec(i)=0.0;
        end
    else
        hs_prec(i)=0.0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %   PRECIPITATION IN BUCKET MODEL
    
    if hs_prec_bucket(i-1)>0.01
        delta_bucket(i)=hs_prec_bucket(i-1)- ZERO;
        hs_prec_bucket(i)=hs_prec(i)+ ZERO;
    elseif hs_prec_bucket(i-1)<=0.01
        delta_bucket(i)=0.0;
        hs_prec_bucket(i)=hs_prec_bucket(i-1)+hs_prec(i);
    end
    
    ros_av(i)=(hs(i-1).*ros(i)+hs_prec_bucket(i).*ros_new(i))./(hs(i-1)+hs_prec_bucket(i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	THERMAL CONDUCTIVITY OF NEW/OLD SNOW
    
    ks_new(i)=(ros_new(i)*10^-3).^2*2.85;   %    (after Abel, 1892)
    ks(i)=(ros(i)*10^-3).^2*2.85;               %    (after Abel, 1892)
   
    ks_av(i)=(ks(i)*hs(i-1)+ks_new(i)*hs_prec_bucket(i))./(hs(i-1)+hs_prec_bucket(i));
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %   ALBEDO COMPUTATION	% (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
       
    if (T0(i-1)<Tfr)
        alpha_i=max(alpha_ow, c11*hi(i-1)^0.28 +0.08);
        alpha_s=0.75;
        alpha_mi=0.70;                                      % (after Perovich, 1996 for compacted snow)
    else
        alpha_i=min(alpha_mi,c12*hi(i-1).^2 +alpha_ow);
        %alpha_s=0.75;                                      % (after Pirazzini et al., 2006 for landfast sea ice of the Baltic)
        alpha_s=0.65;      
        alpha_mi=0.56;                                      % (after Perovich, 1996 for melting white ice)
    end
    
    if (hs(i-1)+hs_prec_bucket(i))>0.1
        alpha(i)=alpha_s;
    elseif hmi(i-1)>hmi_min
        alpha(i)=alpha_mi;
    elseif hi(i-1)>hi_min
        alpha(i)=min(alpha_s,alpha_i+(hs(i-1)+hs_prec_bucket(i)).*(alpha_s-alpha_i)/c10);
    else
        alpha(i)=alpha_ow;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   EXTICTION COEFFICIENT COMPUTATION

    if (T0(i-1)<Tfr)
        ks_daily_severeear(i)=25;
        ks_cloudy(i)=25;
        ks_snow(i)=25;

        kmi_clear(i)=21.05;
        kmi_cloudy(i)=17.75;
        kmi_av(i)=-3.3*Cl(i)+21.05;

        ksi_10_clear(i)=17.1;
        ksi_10_cloudy(i)=10.5;
        ksi_10_av(i)=-6.6*Cl(i)+17.1;

        ksi_clear(i)=1.6;
        ksi_cloudy(i)=1.5;
        ksi_av(i)=-0.1*Cl(i)+1.6;

    elseif (T0(i-1)>=Tfr)
        ks_clear(i)=15;
        ks_cloudy(i)=15;
        ks_snow(i)=15;

        kmi_clear(i)=11.7;
        kmi_cloudy(i)=9.8;
        kmi_av(i)=-1.9*Cl(i)+11.7;

        ksi_10_clear(i)=8.4;
        ksi_10_cloudy(i)=4.6;
        ksi_10_av(i)=-3.8*Cl(i)+8.4;

        ksi_clear(i)=1.5;
        ksi_cloudy(i)=1.4;
        ksi_av(i)=-0.1*Cl(i)+1.5;

    else
        'ERROR EXTICTION COEFFICIENT'
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 	HEAT FLUXES NOT DEPENDENT ON T0

    % 	SHORT WAVE
    Fs(i)=(1 - alpha(i)).*Fsd_cloud(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   LATENT HEAT
    ea(i)=qa(i)*P0/emp;
    
    es(i)=qs(i)*P0/emp;
    Fla(i)=cail*roa*Lai*Ua(i).*(qs(i)-qa(i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%DOWNWARD LONG WAVE
    if hi(i-1)>hi_min
        Fld(i)=((sigma*(Ta(i).^4))-85.6).*(1 + 0.26*Cl(i));                                         % (after Guest, 1997 for Antartica)
    else  
        Fld(i)= (sigma*Ta(i).^4*(0.653 - 0.00535*ea(i))).*(1 + 0.1762*Cl(i).^2);                    % (after Bignami et al, 1995 for Mediterranean Sea)
        Fld(i)=(sigma*Ta(i).^4*(0.685 + 0.00452*ea(i))).*(1 + 0.36*Cl(i).^2);                      % (after Zapadka et al, 2007 for the Baltic)
   end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	OCEANIC FLUXES
    Fw(i)=Fwater;
    %   Fw(i)=row*cpw*cwih*0.05*0.2;%(Tw(i-1)-Tfr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h(i)=h_mix-hi(i-1);

    if hi(i-1)>hi_min
        ice_fr(i)=1.0;
    else
        ice_fr(i)=0.0;
    end

    hs_tot(i)=hs(i-1)+hs_prec_bucket(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %	SURFACE TEMPERATURE ITERATION AND FLUXES

    T0_star(1)=T0(i-1);
    Ts_star(1)=Ts(i-1);
    Tsi_star(1)=Tsi(i-1);
    Ti_star(1)=Ti(i-1);

    kice_surf(i)= k0 + +beta2*Si(i-1)./(Ti(i-1)-273.15);                                                            %   thermal conductivity at the upper boundary of sea ice
    kice_bott(i)=k0 +beta2*Sice_bott(i-1)/(Tfr-273.15);                                                             %   thermal conductivity at the lower boundary of sea ice 
    kice_bio(i)=k0+beta2*Sice_5(i-1)/(Tice_5(i-1));
    ksnow_interm(i)=(kmi*ks_av(i).*(hs(i-1)+hmi(i-1)))./(hs(i-1).*kmi + hmi(i-1).*ks_av(i));                        %   thermalconductivity at the snow/interm layer interface 
    kinterm_ice(i)=(kice_surf(i).*kmi*(hmi(i-1)+hi(i-1)))./(hmi(i-1).*kice_surf(i)+hi(i-1).*kmi);                   %   thermalconductivity at the interm layer/sea ice interface 
    ksnow_ice(i)=(kice_surf(i).*ks_av(i).*(hs(i-1)+hi(i-1)))./(hs(i-1).*kice_surf(i) + hi(i-1).*ks_av(i));          %   thermal conductivity at the snow/sea ice layer interface 
  
    Kice_bott(i)=2*kice_bott(i)./hi(i-1);
    Kice_surf(i)=2*kice_surf(i)./hi(i-1);
    Kice_bio(i)=2*kice_bio(i)./hi_bio(i-1);
    Kbott_bio(i)=2*kice_bott(i)./hi_bio(i-1);
    
    Kmi(i)=2*kmi/(hmi(i-1));
    Ks(i)=2*ks_av(i)./hs(i-1);
    Ksnow_interm(i)=2*ksnow_interm(i)./(hmi(i-1)+hs(i-1));
    Kinterm_ice(i)=2*kinterm_ice(i)./(hmi(i-1)+hi(i-1));
    Ksnow_ice(i)=2*ksnow_ice(i)./(hi(i-1)+hs(i-1));
    
    %mus(i)=deltat/(c0*(ros(i).*hs(i-1)));
    mus(i)=deltat/(c0*(ros_new(i)*hs_prec_bucket(i-1)+ros(i).*hs(i-1)));
    mumi(i)=deltat/(romi*c0*hmi(i-1));
    cpi_bio(i)=c0 + (L0s *mu*Sice_bio(i-1))/(Tice_bio(i-1)^2);
    mui_bio(i)=deltat/(ro_sice_bulk(i-1)*10^3.*cpi_bio(i).*hi_bio(i-1));
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%% ITERATION START %%%%%%%%%%%%%%%%%%%%
    for j=2:100

        Flu_it=emis*sigma*T0_star(j-1).^4;                                                      %   UPWARD LONG WAVE
        Fse_it=roa*cpair*cais*Ua(i).*(T0_star(j-1)-Ta(i));                                      %   SEMSIBLE HEAT
        Fl_it= Flu_it-Fld(i);                                                                   %   NET LONG WAVE
        F_it= - Fs(i) + Fl_it + Fla(i) + Fse_it;                                                %   NET SURFACE FLUXES

        la=4*sigma*T0_star(j-1).^3;                                                             %   LONG WAVE (differential) upward ==>negative downward
        sen=roa*cpair*cais*Ua(i);                                                               %   SENSIBLE HEAT (differential) upward ==>negative downward
        lat=cail*roa*Lai*Ua(i).*((c1*(Tb-c2)*qs(i))./((T0_star(j-1)+Tb-c2).^2));                %   LATENT HEAT (differential) upward ==>negative downward
        cons=ks_av(i)./(hs(i-1)/2);                                                             %   SNOW CONDUCTIVE (differential)
        coni=kice_surf(i)./(hi(i-1)/2);                                                         %   SEA ICE CONDUCTIVE (differential)
        conmi=kmi/(hmi(i-1)/2);                                                                 %   SNOW ICE/SUPERIMPOSED ICE (differential)
        cpi_it=c0 + ((L0s *mu*Si(i-1))/((Ti_star(j-1)-273.15)*(Ti_star(1)-273.15)));
        mui=deltat/(ro_sice_bulk(i-1)*10^3.*cpi_it.*hi(i-1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %    TEMPERATURE
        if ((hs(i-1)+hs_prec_bucket(i))>hs_min)
             
            %   SEA ICE MODEL ON/SLAB OCEAN OFF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if hmi(i-1)>hmi_min

                %   CASE 1: 3 LAYERS (SNOW + INTERM LAYER + SEA ICE)

                I0_it=vis_fr*Fs(i).*expm(-ks_snow(i).*(hs(i-1)+hs_prec_bucket(i)));
                IM_it=I0_it.*expm(-kmi_av(i).*hmi(i-1));
                if hi(i-1)>=0.1
                    ISI_10_it=IM_it.*expm(-ksi_10_av(i)*0.1);
                    ISI_it=ISI_10_it.*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10_it=IM_it.*expm(-ksi_10_av(i).*hi(i-1));
                    ISI_it=ISI_10_it;
                end

                A=[(1-Ks(i)./(la+sen+lat+cons)) (Ks(i)./(la+sen+lat+cons)) 0 0;
                    (mus(i).*Ks(i)) (1-mus(i).*(Ksnow_interm(i)+Ks(i))) (mus(i).*Ksnow_interm(i)) 0;
                    0 (mumi(i).*Ksnow_interm(i)) (1- mumi(i).*(Ksnow_interm(i)+Kinterm_ice(i))) (mumi(i)*Kinterm_ice(i));
                    0 0 (mui*Kinterm_ice(i)) (1-mui.*(Kinterm_ice(i)+Kice_bott(i)))];

                b=[((-F_it-I0_it)/(la+sen+lat+cons)) (I0_it*mus(i)) (IM_it*mumi(i)) (mui.*Kice_bott(i)*Tfr + mui.*ISI_it)]';
                T_prev=[T0_star(j-1) Ts_star(1) Tsi_star(1) Ti_star(1)]';
                T(:,j)=A*T_prev + b;
                T0_star(j)=T(1,j);
                Ts_star(j)=T(2,j);
                Tsi_star(j)=T(3,j);
                Ti_star(j)=T(4,j);

                T0_iter=T0_star(end);
                Ts_iter=Ts_star(end);
                Tsi_iter=Tsi_star(end);
                Ti_iter=Ti_star(end);
                Tmix_star=Tfr;

                if Ts_iter>Tfr
                    Ts_iter=Tfr;
                end
                
                if Tsi_iter>Tfr
                     Tsi_iter=Tfr;
                     Ts_iter=Tfr;
                end
                if Ti_iter>Tfr
                    Ti_iter=Tfr;
                    Tsi_iter=Tfr;
                    Ts_iter=Tfr;
                end

                if T0_iter>Tfr
                    T0_iter=Tfr;
                    B=[0 0 0 0;
                        (mus(i).*Ks(i)) (1-mus(i).*(Ksnow_interm(i)+Ks(i))) (mus(i).*Ksnow_interm(i)) 0;
                        0 (mumi(i).*Ksnow_interm(i)) (1- mumi(i).*(Ksnow_interm(i)+Kinterm_ice(i))) (mumi(i).*Kinterm_ice(i));
                        0 0 (mui.*Kinterm_ice(i)) (1-mui.*(Kinterm_ice(i)+Kice_bott(i)))];

                    c=[0 (I0_it*mus(i)) (IM_it*mumi(i)) (mui.*Kice_bott(i)*Tfr + mui.*ISI_it)]';
                    T_prev_melt=[Tfr Ts(i-1) Tsi(i-1) Ti(i-1)]';
                    T_melt=B*T_prev_melt + c;
                    Ts_iter=T_melt(2);
                    if Ts_iter<Tfr
                        Tsi_iter=T_melt(3);
                        if Tsi_iter<Tfr
                            Ti_iter=T_melt(4);
                            if Ti_iter>Tfr
                                Ti_iter=Tfr;
                            end
                        else
                            Tsi_iter=Tfr;
                            Ti_iter=(mui.*Kinterm_ice(i)).*Tfr + (1-mui.*(Kinterm_ice(i)+Kice_bott(i))).*Ti(i-1) + (mui.*Kice_bott(i).*Tfr);
                            if Ti_iter>Tfr
                                Ti_iter=Tfr;
                            end
                        end
                    else
                        Ts_iter=Tfr;
                        C=[0 0 0 0;
                            0 0 0 0;
                            0 (mumi(i).*Ksnow_interm(i)) (1- mumi(i).*(Ksnow_interm(i)+Kinterm_ice(i))) (mumi(i).*Kinterm_ice(i));
                            0 0 (mui.*Kinterm_ice(i)) (1-mui.*(Kinterm_ice(i)+Kice_bott(i)))];

                        d=[0 0 (IM_it*mumi(i)) (mui.*Kice_bott(i)*Tfr + mui*ISI_it)]';
                        T_prev_melt_2=[Tfr Tfr Tsi(i-1) Ti(i-1)]';
                        T_melt_2=C*T_prev_melt_2 + d;
                        Tsi_iter=T_melt_2(3);
                        if Tsi_iter<Tfr
                            Ti_iter=T_melt_2(4);
                            if Ti_iter>Tfr
                                Ti_iter=Tfr;
                            end
                        else
                            Tsi_iter=Tfr;
                            Ti_iter=(mui.*Kinterm_ice(i)).*Tfr + (1-mui.*(Kinterm_ice(i)+Kice_bott(i))).*Ti(i-1) + (mui.*Kice_bott(i).*Tfr + mui.*ISI_it);
                            if Ti_iter>Tfr
                                Ti_iter=Tfr;
                            end
                        end
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif hmi(i-1)<=hmi_min

                %   CASE 2: 2 LAYERS (SNOW + SEA ICE)

                I0_it=vis_fr*Fs(i).*expm(-ks_snow(i).*(hs(i-1)+hs_prec_bucket(i)));
                IM_it=I0_it;
                if hi(i-1)>=0.1
                    ISI_10_it=IM_it.*expm(-ksi_10_av(i)*0.1);
                    ISI_it=ISI_10_it.*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10_it=IM_it.*expm(-ksi_10_av(i).*hi(i-1));
                    ISI_it=ISI_10_it;
                end

                A=[(1-Ks(i)./(la+sen+lat+cons)) (Ks(i)./(la+sen+lat+cons)) 0;
                    (mus(i).*Ks(i)) (1- mus(i).*(Ksnow_ice(i)+Ks(i))) (mus(i).*Ksnow_ice(i));
                    0 (mui.*Ksnow_ice(i)) (1 - mui.*(Ksnow_ice(i) + Kice_bott(i)))];

                b=[((-F_it-I0_it)/(la+sen+lat+cons)) (I0_it*mus(i)) (mui.*Kice_bott(i).*Tfr + mui.*ISI_it)]';
                T_prev=[T0_star(j-1) Ts_star(1) Ti_star(1)]';
                T(:,j)=A*T_prev + b;
                T0_star(j)=T(1,j);
                Ts_star(j)=T(2,j);
                Ti_star(j)=T(3,j);

                T0_iter=T0_star(end);
                Ts_iter=Ts_star(end);
                Ti_iter=Ti_star(end);
                Tmix_star=Tfr;

                if Ts_iter>Tfr
                    Ts_iter=Tfr;
                end
                if Ti_iter>Tfr
                    Ti_iter=Tfr;
                    Ts_iter=Tfr;
                end
                

                if T0_iter>Tfr
                    T0_iter=Tfr;
                    B=[0 0 0;
                        (mus(i).*Ks(i)) (1- mus(i).*(Ksnow_ice(i)+Ks(i))) (mus(i).*Ksnow_ice(i));
                        0 (mui.*Ksnow_ice(i)) (1 - mui.*(Ksnow_ice(i) + Kice_bott(i)))];

                    c=[0 (I0_it*mus(i)) (mui.*Kice_bott(i).*Tfr + mui.*ISI_it)]';
                    T_prev_melt=[Tfr Ts(i-1) Ti(i-1)]';
                    T_melt=B*T_prev_melt + c;
                    Ts_iter=T_melt(2);
                    if Ts_iter<Tfr
                        Ti_iter=T_melt(3);
                    else
                        Ts_iter=Tfr;
                        Ti_iter=(mui.*Ksnow_ice(i))*Tfr + (1 - mui.*(Ksnow_ice(i) + Kice_bott(i)))*Ti(i-1) + (mui.*Kice_bott(i).*Tfr + mui.*ISI_it);
                        if Ti_iter>Tfr
                            Ti_iter=Tfr;
                        end
                    end
                end
                    Tsi_iter=Tfr;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  ((hs(i-1)+hs_prec_bucket(i))<=hs_min)

            if hmi(i-1)>hmi_min

                %   CASE 3: 2 LAYERS (INTERM LAYER + SEA ICE)

                I0_it=vis_fr*Fs(i).*expm(-kmi_av(i).*hmi(i-1));
                IM_it=I0_it;
                if hi(i-1)>=0.1
                    ISI_10_it=IM_it.*expm(-ksi_10_av(i)*0.1);
                    ISI_it=ISI_10_it.*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10_it=IM_it.*expm(-ksi_10_av(i).*hi(i-1));
                    ISI_it=ISI_10_it;
                end

                A=[(1-Kmi(i)./(la+sen+lat+conmi)) (Kmi(i)./(la+sen+lat+conmi)) 0;
                    (mumi(i).*Kmi(i)) (1 - mumi(i).*(Kinterm_ice(i) +Kmi(i))) (mumi(i).*Kinterm_ice(i));
                    0 (mui.*Kinterm_ice(i)) (1 - mui.*(Kinterm_ice(i) + Kice_bott(i)))];

                b=[((-F_it-I0_it)/(la+sen+lat+conmi)) (I0_it*mumi(i)) (mui.*Kice_bott(i).*Tfr + mui.*ISI_it)]';
                T_prev=[T0_star(j-1) Tsi_star(1) Ti_star(1)]';
                T(:,j)=A*T_prev + b;
                T0_star(j)=T(1,j);
                Tsi_star(j)=T(2,j);
                Ti_star(j)=T(3,j);
                
                %Ti_bio_prova(i)=Ti_bio_prova(i-1)*((mui_bio(i)*Kice(i)+()(Io-bio(i)*mui_bio(i)

                T0_iter=T0_star(end);
                Tsi_iter=Tsi_star(end);
                Ti_iter=Ti_star(end);
                Tmix_star=Tfr;

                if Tsi_iter>Tfr
                    Tsi_iter=Tfr;
                end
                if Ti_iter>Tfr
                    Ti_iter=Tfr;
                    Tsi_iter=Tfr;
                end

                if T0_iter>Tfr
                    T0_iter=Tfr;
                    B=[(1-Kmi(i)./(la+sen+lat+conmi)) (Kmi(i)./(la+sen+lat+conmi)) 0;
                        (mumi(i).*Kmi(i)) (1 - mumi(i).*(Kinterm_ice(i) +Kmi(i))) (mumi(i).*Kinterm_ice(i));
                        0 (mui.*Kinterm_ice(i)) (1 - mui.*(Kinterm_ice(i) + Kice_bott(i)))];

                    c=[0 (I0_it*mumi(i)) (mui.*Kice_bott(i).*Tfr + mui.*ISI_it)]';
                    T_prev_melt=[Tfr Tsi(i-1) Ti(i-1)]';
                    T_melt=B*T_prev_melt + c;
                    Tsi_iter=T_melt(2);
                    if Tsi_iter<Tfr
                        Ti_iter=T_melt(3);
                    else
                        Tsi_iter=Tfr;
                        Ti_iter=(mui.*Kinterm_ice(i)).*Tfr + (1 - mui.*(Kinterm_ice(i) + Kice_bott(i))).*Ti(i-1) + (mui.*Kice_bott(i).*Tfr + mui.*ISI_it);
                        if Ti>Tfr
                            Ti=Tfr;
                        end
                    end
                end             
                    Ts_iter=Tfr;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif hmi(i-1)<=hmi_min
                if hi(i-1)>hi_min
                    %   CASE 4: ONLY SEA ICE

                    if hi(i-1)>=0.1
                        ISI_10_it=vis_fr*Fs(i).*expm(-ksi_10_av(i)*0.1);
                        ISI_it=ISI_10_it.*expm(-ksi_av(i).*(hi(i-1)-0.1));
                    elseif hi(i-1)<0.1
                        ISI_10_it=vis_fr*Fs(i).*expm(-ksi_10_av(i).*hi(i-1));
                        ISI_it=ISI_10_it;
                    end

                    I0_it=ISI_10_it;
                    IM_it=ISI_10_it;

                    A=[(1-Kice_surf(i)./(la+sen+lat+coni)) (Kice_surf(i)./(la+sen+lat+coni)) 0;
                        (mui.*Kice_surf(i)) (1 - mui.*(Kice_surf(i) + Kice_bott(i))) 0;
                        0 0 0];

                    b=[((-F_it-I0_it)/(la+sen+lat+coni)) (mui.*Kice_bott(i).*Tfr + mui.*ISI_it) 0]';
                    T_prev=[T0_star(j-1) Ti_star(1) Ti_star(1)]';
                    T(:,j)=A*T_prev + b;
                    T0_star(j)=T(1,j);
                    Ti_star(j)=T(2,j);

                    T0_iter=T0_star(j);
                    Ti_iter=Ti_star(j);
                    Tmix_star=Tfr;

                    if Ti_iter>Tfr
                        Ti_iter=Tfr;
                    end

                    if T0_iter>Tfr
                        T0_iter=Tfr;
                        Ti_iter=(mui.*Kice_surf(i)).*Tfr + (1 - mui.*(Kice_surf(i) + Kice_bott(i))).*Ti(i-1) + (mui.*Kice_bott(i).*Tfr + mui.*ISI_it);
                        if Ti_iter>Tfr
                            Ti_iter=Tfr;
                        end
                    end
                        Ts_iter=Tfr;
                        Tsi_iter=Tfr;
                                            
                elseif hi(i-1)<=hi_min
                    
                    %   SEA ICE MODEL OFF/SLAB OCEAN ON 
                    %   CASE 5: 0 LAYER OR FRAZIL ICE FORMATION
                    
                    I0_it=Fs(i).*((infra_fr*expm(-k_ocean_red*h_mix)) + ((1-infra_fr)*expm(-k_ocean_vis*h_mix)));
                    Tmix_star=Tmix(i-1)-((F_it + I0_it + R(i-1))./(h_mix*row*cpw)).*deltat;
                    T0_star(j)=Tmix_star;
                    Ts_star(j)=Tfr;
                    Tsi_star(j)=Tfr;
                    Ti_star(j)=Tfr;
                    T0_iter=T0_star(j);
                    Ts_iter=Ts_star(j);
                    Tsi_iter=Tsi_star(j);
                    Ti_iter=Ti_star(j);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            else
                'ERROR ITERATION'
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if abs(T0_star(j)-T0_star(j-1))< 0.01, break                    %   CONVERGENCE CRITERION: 0.01 K
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length (T0_star)>20
        'WARNING: ITERATION >20'
    end
    if length (T0_star)==100                                             %   MAX NUMBER OF ITERATIONS FOR CONVERGENCE
        'NOT CONVERGENT'
    end

    if max(size(T0_star))==100, return                                   %   STOP THE SEA ICE MODEL IF NOT CONVERGENT
    end

    T0(i)=T0_iter;                                                       %   NEW SURFACE TEMPERATURE
    Ts(i)=Ts_iter;                                                       %   NEW SNOW TEMPERATURE
    Tsi(i)=Tsi_iter;                                                     %   NEW SNOW ICE/SUPERIMPOSED ICE TEMPERATURE
    Ti(i)=Ti_iter;                                                       %   NEW SEA ICE TEMPERATURE
    Tmix(i)=Tmix_star;                                                   %   NEW MIX LAYER TEMPERATURE
    
    
    if hmi(i-1)>hmi_min
        Tice(i)=((kice_surf(i).*Ti(i)/(hi(i-1)/2) + kmi.*Tsi(i)/(hmi(i-1)/2))./(kice_surf(i)./(hi(i-1)/2) +kmi./(hmi(i-1)/2)))-273.15;
    elseif hmi(i-1)<=hmi_min
        if ((hs(i-1)+hs_prec_bucket(i))>hs_min)
        Tice(i)=((kice_surf(i).*Ti(i)/(hi(i-1)/2) + ks_av(i).*Ts(i)/(hs(i-1)/2))./(kice_surf(i)./(hi(i-1)/2) + ks_av(i)./(hs(i-1)/2)))-273.15;
        elseif ((hs(i-1)+hs_prec_bucket(i))<=hs_min)
            if hi(i-1)>hi_min
                Tice(i)=T0(i)-273.1499;
            elseif hi(i-1)<=hi_min
                Tice(i)=Tfr-273.15;
            end
        end
    end
   
    clear T0_star Ts_star Tsi_star Ti_star Tmix_star T_prev T A a B b C d D e Flu_it Fse_it Fl_it F_it la sen lat ISI_it ISI_10_it IMI_it I0_it
    clear T_prev_melt T_melt T_prev_melt_2 T_melt_2 C c D T0_iter Ts_iter Tsi_iter Ti_iter cons coni conmi mui cpi_it
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   RICOMPUTE SURFACE FLUXES
    
    Flu(i)=emis*sigma*T0(i).^4;                                         %   UPWARD LONGWAVE RADIATION
    Fl(i)= Flu(i)-Fld(i);                                               %   NET LONGWAVE RADIATION
    Fse(i)=roa*cpair*cais*Ua(i).*(T0(i)-Ta(i));                         %   SENSIBLE HE
    F(i)=-Fs(i) + Fl(i) +Fla(i) +Fse(i);                                %   NET SURFACE FLUXES
    R(i)=0.0;                                                           %   INITIAL EXTRA HEAT IS SET TO ZERO
    
    Qs(i)=ros(i).*((c0*(Tfrs-Ts(i))) + L0s);
    Qmi(i)=romi.*((c0*(Tfrs-Tsi(i))) + L0s);
    Qi_surf(i)=ro_sice_surf(i-1).*10^3*((c0*(Tfr-Ti(i))) + L0s*(1+mu*Si(i-1)/(Ti(i)-273.15)));
    Qi_bott(i)=ro_sice_bott(i-1).*10^3*((c0*(Tfr-Ti(i))) + L0s*(1+mu*Si(i-1)/(Ti(i)-273.15)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEA ICE GROWTH/DECAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GROWTH SEASON

    if (T0(i)<Tfr)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if hmi(i-1)>hmi_min

            if (hs(i-1)+hs_prec_bucket(i))>hs_min

                %   CASE 1: SNOW + INTERM LAYER + SEA ICE
                
               
                
                Fsurf(i)=Ks(i).*(Ts(i)-T0(i));
                Fbott(i)=Kice_bott(i)*(Tfr-Ti(i));
                
                I0(i)=vis_fr*Fs(i).*expm(-ks_snow(i).*(hs(i-1)+hs_prec_bucket(i)));
                IM(i)=I0(i).*expm(-kmi_av(i).*hmi(i-1));
                
                if hi(i-1)>=0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i)*0.1);
                    ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i).*hi(i-1));
                    ISI(i)=ISI_10(i);
                end               
                
                delta_sublim(i)=(Fla(i)*deltat)/(ros_new(i).*Lv-Qs(i));                
                deltahs_melt_surf(i)=0.0;
                deltahi_melt_surf(i)=0.0;
                deltahmi_melt_surf(i)=0.0;
                deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);

                if (hs(i-1)+hs_prec_bucket(i)-delta_sublim(i))> ((hi(i-1).*(row-ro_sice_surf(i-1)*10^3)/ros_av(i)) + (hmi(i-1)*(row-romi)/ros_av(i))) %   SNOW ICE FORMATION
                    hs(i)=hs(i-1)+delta_bucket(i) - beta*((-row*hi(i-1) + ros(i).*(hs(i-1)-delta_sublim(i))+ 10^3*ro_sice_surf(i-1).*hi(i-1) - row*hmi(i-1) + romi*hmi(i-1) + ros_new(i)*hs_prec_bucket(i))/(-romi + row+ beta*ros_av(i)));
                    hmi_new(i)=((-row*hi(i-1) + ros(i).*(hs(i-1)-delta_sublim(i)) + 10^3*ro_sice_surf(i-1).*hi(i-1) - row*hmi(i-1) + romi*hmi(i-1)+ ros_new(i)*hs_prec_bucket(i))/(-romi + row+ beta*ros_av(i))*si_fract);
                    hmi(i)=hmi(i-1) + hmi_new(i);
                else
                    hs(i)=hs(i-1)+ delta_bucket(i) - delta_sublim(i);
                    hmi_new(i)=0.0;
                    hmi(i)=hmi(i-1);
                end
                hi(i)=hi(i-1)+ deltahi_bott(i);
                
                if hs(i)<hs_min
                    hs(i)=hs_prec_bucket(i);
                    hs_prec_bucket(i)=0.0;                   
                end
                    
            elseif (hs(i-1)+hs_prec_bucket(i))<=hs_min

                %   CASE 2: INTERM LAYER + SEA ICE

                Fsurf(i)=Kmi(i).*(Tsi(i)-T0(i));
                Fbott(i)=Kice_bott(i).*(Tfr-Ti(i));
                I0(i)=vis_fr*Fs(i).*expm(-kmi_av(i).*hmi(i-1));
                IM(i)=I0(i);

                if hi(i-1)>=0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i)*0.1);
                    ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i).*hi(i-1));
                    ISI(i)=ISI_10(i);
                end
                
                delta_sublim(i)=(Fla(i)*deltat)/(romi.*Lv-Qmi(i));
                deltahmi_melt_surf(i)=0.0;
                deltahi_melt_surf(i)=0.0;
                deltahs_melt_surf(i)=0.0;
                deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);

                hs(i)=hs(i-1)+delta_bucket(i);
                hmi_new(i)=0.0;
                hmi(i)=hmi(i-1)-delta_sublim(i);
                hi(i)=hi(i-1)+ deltahi_bott(i);
                %fb(i)=0.0;
                %hw(i)=0.0;
                %Ssnowice(i)=Ssnowice(i-1);
            end
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        elseif hmi(i-1)<=hmi_min

            if (hs(i-1)+hs_prec_bucket(i))>hs_min
                
                %   CASE 3: SNOW + SEA ICE
                
                Fsurf(i)=Ks(i).*(Ts(i)-T0(i));
                Fbott(i)=Kice_bott(i).*(Tfr-Ti(i));

                I0(i)=vis_fr*Fs(i).*expm(-ks_snow(i).*(hs(i-1)+hs_prec_bucket(i)));
                IM(i)=I0(i);

                if hi(i-1)>=0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i)*0.1);
                    ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i).*hi(i-1));
                    ISI(i)=ISI_10(i);
                end
                
                delta_sublim(i)=(Fla(i)*deltat)/(ros_new(i).*Lv-Qs(i));                
                deltahs_melt_surf(i)=0.0;
                deltahi_melt_surf(i)=0.0;
                deltahmi_melt_surf(i)=0.0;
                deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);

                if (hs(i-1) + hs_prec_bucket(i) - delta_sublim(i))> (hi(i-1).*(row-ro_sice_surf(i-1)*10^3)/ros_av(i) - hs_prec_bucket(i).*ros_new(i)/ros(i))  %   SNOW ICE FORMATION
                    hs(i)=hs(i-1) + delta_bucket(i) -delta_sublim(i) -beta*((-row*hi(i-1) + ros(i).*hs(i-1)+ 10^3*ro_sice_surf(i-1).*hi(i-1) + ros_new(i).*hs_prec_bucket(i))/(-romi + row+ beta*ros_av(i)));
                    hmi_new(i)=((-row*hi(i-1) + ros(i).*(hs(i-1)-delta_sublim(i))+ 10^3*ro_sice_surf(i-1).*hi(i-1) + ros_new(i)*hs_prec_bucket(i))/(-romi + row+ beta*ros_av(i)))*si_fract;
                    hmi(i)=hmi(i-1) + hmi_new(i);                  
                else                   
                    hs(i)=hs(i-1) + delta_bucket(i) - delta_sublim(i);
                    hmi_new(i)=0.0;
                    hmi(i)=hmi(i-1);
                end
                
                if hs(i)<ZERO
                    hs(i)=hs_prec_bucket(i);
                    hs_prec_bucket(i)=0.0;
                end
                
                hi(i)= hi(i-1) + deltahi_bott(i);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            elseif ((hs(i-1) + hs_prec_bucket(i))<=hs_min)
                
                %   CASE 4: ONLY SEA ICE
                Fsurf(i)=Kice_surf(i).*(Ti(i)-T0(i));
                Fbott(i)=Kice_bott(i).*(Tfr-Ti(i));

                if hi(i-1)>=0.1
                    ISI_10(i)=vis_fr*Fs(i).*expm(-ksi_10_av(i)*0.1);
                    ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10(i)=vis_fr*Fs(i).*expm(-ksi_10_av(i).*hi(i-1));
                    ISI(i)=ISI_10(i);
                end

                I0(i)=ISI_10(i);
                IM(i)=ISI_10(i);
                
                if hi(i-1)>hi_min
                    delta_sublim(i)=(Fla(i)*deltat)/(ro_sice_surf(i-1).*Lv-Qi_surf(i));
                    deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);
                else
                    delta_sublim(i)=0.0;
                    ISI(i)=Fs(i).*((infra_fr*expm(-k_ocean_red*h_mix)) + ((1-infra_fr)*expm(-k_ocean_vis*h_mix)));
                    deltahi_bott(i)=deltat*(F(i)+ISI(i))./Qi_bott(i);
                end

                deltahs_melt_surf(i)=0.0;
                deltahmi_melt_surf(i)=0.0;
                deltahi_melt_surf(i)=0.0; 

                hs(i)=hs(i-1)+delta_bucket(i);
                hmi_new(i)=0.0;
                hmi(i)=hmi(i-1);
                hi(i)=  hi(i-1) + deltahi_bott(i) - delta_sublim(i);
                
                if hi(i)<ZERO
                    R(i)=(ro_sice_bott(i-1).*10^3*((c0*(Tfr-Ti(i))) + L0s*(1+Si(i-1)/(Ti(i)-273.15))))*(hi(i)-ZERO)./deltat;  % EXTRA HEAT FOR ENERGY CONSERVATION
                    hi(i)=ZERO;
                    hs(i)=ZERO;
                    hmi(i)=ZERO;
                end
            end
        else
            'ERROR GROWTH SEASON'
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HALODYNAMIC SUBMODEL FOR THE GROWTH SEASON
        if hi(i)>ZERO

            vi(i)=deltahi_bott(i)./deltat;
               %keff(i)=0.113/(0.113 + 0.887*expm(-2.66*10^4*vi(i))); % for Baltic (Gransgok)
            keff(i) = 0.12/(0.12+0.88*expm(-4.2*10^4*vi(i))); % with 4.2*10^4 in s/cm for normal salinity sea water >30, it is usuallly used the parameterization by Nawako and Shina, 1983
            Sice_bott(i)=keff(i).*Sw(i);
            Sbr_bott(i)=-(Tfr-273.1499)/mu;
            Vbr_bott(i)=Sice_bott(i)./Sbr_bott(i);
            Q(i)=0.0;
            
             %   SNOW ICE LAYER
             if hmi_new(i)>0.0
                snow_fr(i)=beta/(ro_sice_surf(i-1)*10^3./ros_av(i));
                sea_water_fr(i)=1-snow_fr(i);
                Ssnowice_new(i)= (hmi_new(i).*sea_water_fr(i).*Sw(i))./hmi_new(i);
                Ssnowice(i)= (Ssnowice(i-1)*hmi(i-1)+Ssnowice_new(i)*hmi_new(i))./hmi(i);
             else 
                 snow_fr(i)=snow_fr(i-1);
                 sea_water_fr(i)=sea_water_fr(i-1);
                 Ssnowice_new(i)=0.0;
                 Ssnowice(i)=Ssnowice(i-1);
             end
%                                             
%             Ms(i)=ros(i).*hs(i)+hs_prec_bucket(i).*ros_new(i);         %   snow mass
%             Mmi(i)=romi*hmi(i);                                         %   snow ice mass
%             hw(i)=(Mmi(i)+ Ms(i))/row;                                  %   water equivalent thickness
%             Mw(i)=si_fract*row*hw(i);                                            %   sea water mass      
%             if hmi(i)>hmi_min
%                 Ssnowice(i)=(Ssnow*Ms(i)+Sw(i),*Mw(i))./(Ms(i)+Mw(i));      %   snow ice salinity
%                 %Ssnowice(i)=Ssnowice(i-1) + (Ssnow*(Ms(i)-Ms(i-1))+Sw(i),*(Mw(i)-Mw(i-1))./(Ms(i)+Mw(i)));      %   snow ice salinity
%             else
%                 Ssnowice(i)=0.0;                        
%             end
%                 

            if (Vbr_ice(i-1)> 0.05)                                                         %   DESALINATION
                Sice_mix(i)=((hi(i)-hi(i-1)).*Sice_bott(i) + hi(i-1).*Sice(i-1))./hi(i);
                Sice(i) = Sice_mix(i) + (rosa*(1-viola*Vbr_ice(i-1)).*(Tfr-273.15-Tice(i)))*deltat;
                Sbr_ice(i)=-Tice(i)/mu;
                Vbr_ice(i)=Sice(i)./Sbr_ice(i);

                Sice_5(i)=Sice(i);                                                          %   DESALINATION=0, IT STARTS HERE
                Sbr_5(i)=Sbr_ice(i);
                Tice_5(i)=Tice(i);
%                 hi_5(i)=hi(i);
                
                Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.1499) + ki_ice_bio(i-1)*Tice(i))/(ki_bio_bott(i-1)+ki_ice_bott(i-1));
                
            elseif (Vbr_ice(i-1)<=0.05)                                                    %   NO DESALINATION
                Sice_mix(i)=((hi(i)-hi(i-1)).*Sice_bott(i) + hi(i-1).*Sice(i-1))/hi(i);
                Sice(i)=Sice_mix(i);
                Sbr_ice(i)=-Tice(i)/mu;
                Vbr_ice(i)=Sice(i)./Sbr_ice(i);    
                
                
                Sice_5_mix(i)=((hi_5(i-1)*Sice_5(i-1)+(hi(i)-hi(i-1)).*Sice_bott(i)))/(hi(i)-hi(i-1)+hi_5(i-1));
                Sice_5(i)=Sice_5_mix(i) + (rosa*(1-viola*0.0499)*(Tfr-273.15-Tice_5(i-1)))*deltat;
                Sbr_5(i)=Sice_5(i)/0.0499;
                Tice_5(i)=-mu*Sbr_5(i);
%                 if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))               
%                     hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.1499)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.1499))));
%                 else
%                     %hi_5(i)=ZERO;
%                     hi_5(i)=hi_5(i-1);
%                 end
%                 
                Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.1499) + ki_ice_5(i-1)*Tice_5(i))/(ki_bio_bott(i-1)+ki_ice_5(i-1));
                
            end
            
            hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.1499)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.1499))));              
            hi_bio(i)=hi_5(i);           
            Sice_bio_mix(i)=((hi(i)-hi(i-1))*Sice_bott(i) + hi_bio(i-1)*Sice_bio(i-1))./(hi(i)-hi(i-1)+hi_bio(i-1));  % DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
            Sice_bio(i)=Sice_bio_mix(i)+(rosa*(1-viola*Vbr_bio(i-1)).*(Tfr-273.15-Tice_bio(i)))*deltat;
            Sbr_bio(i)=-Tice_bio(i)/mu;
            Vbr_bio(i)=Sice_bio(i)./Sbr_bio(i);
            
            Si_mix(i)=((hi(i)-hi(i-1)).*Sice_bott(i) + (hi(i-1)/2).*Si(i-1))./(hi(i)-hi(i-1)+hi(i-1)/2);
            if Vbr_i(i-1)<=0.05               
                Si(i)=Si_mix(i);               
            else
                Si(i)=Si_mix(i) + (rosa*(1-viola*Vbr_i(i-1)).*(273.1499-Ti(i)))*deltat;
            end
            Sbr_i(i)=-(Ti(i)-273.1499)/mu;
            Vbr_i(i)=Si(i)./Sbr_i(i);
            
            
        elseif hi(i)<=ZERO
            % NO SALINITY COMPUTATION

            snow_fr(i)=0.0;
            sea_water_fr(i)=0.0;
            Ssnowice_new(i)=0.0;
            Ssnowice(i)=0.0; 
            
            Q(i)=ZERO;
            Sice(i)=Sw(i);
            Sice_bott(i)=Sw(i);
            
            Sice_mix(i)=Sw(i);
            Sice_5_mix(i)=Sw(i);
            Sice_bio_mix(i)=Sw(i);
        
            Sbr_ice(i)=Sw(i);
            Vbr_ice(i)=1.0;

            Sbr_5(i)=Sw(i);
            Sice_5(i)=Sw(i);
            Tice_5(i)=Tfr-273.15;
            hi_5(i)=ZERO;
            
            Si(i)=Sw(i);
            Vbr_i(i)=1.0;
            Sbr_i(i)=Sw(i);
            
            Sbr_bio(i)=Sw(i);
            Sice_bio(i)=Sw(i);
            Tice_bio(i)=Tfr-273.15;
            Vbr_bio(i)=1.0;
            hi_bio(i)=hi_5(i);
            ISI_bio(i)=0.0;
            
            Sice_bott(i)=Sw(i);
            Sbr_bott(i)=Sw(i);
            Vbr_bott(i)=1.0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     END OF THE GROWING SEASON ---- START THE MELTING SEASON     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %   MELT SEASON

        %   THERMODYNAMIC AND HALODYNAMIC SUBMODEL

    elseif (T0(i)>=Tfr)

        keff(i)=1.000;
        vi(i)=0.0;        
        hmi_new(i)=0.0;
        snow_fr(i)=0.0;
        sea_water_fr(i)=0.0;
        Ssnowice_new(i)=0.0;
        
        if ((hs(i-1)+ hs_prec_bucket(i))>hs_min)
            
            %    CASE 1: MELTING SNOW
            T0(i)=Tfr;
            I0(i)=vis_fr*Fs(i).*expm(-ks_snow(i).*(hs(i-1)+hs_prec_bucket(i)));
            Fsurf(i)=Ks(i).*(Ts(i)-T0(i));
            Fbott(i)=Kice_bott(i).*(Tfr-Ti(i));
            
            if hmi(i-1)>hmi_min                
                IM(i)=I0(i).*expm(-kmi_av(i).*hmi(i-1));
            elseif hmi(i-1)<=hmi_min  
                hmi(i)=ZERO;
                IM(i)=I0(i);
            end

            if hi(i-1)>=0.1
                ISI_10(i)=IM(i).*expm(-ksi_10_av(i)*0.1);
                ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
            elseif hi(i-1)<0.1
                ISI_10(i)=IM(i).*expm(-ksi_10_av(i).*hi(i-1));
                ISI(i)=ISI_10(i);
            end

            
            if Fsurf(i)>=(F(i)+I0(i))
                deltahs_melt_surf(i)=deltat*(F(i)+I0(i)-Fsurf(i))./Qs(i);
            else
                deltahs_melt_surf(i)=0.0;
            end
            
            delta_sublim(i)=(Fla(i)*deltat)/(ros_av(i).*Lv-Qs(i));
            deltahi_melt_surf(i)=0.0;
            deltahmi_melt_surf(i)=0.0;
            deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);

            hs(i)=hs(i-1)+deltahs_melt_surf(i)+delta_bucket(i)-delta_sublim(i);
            hmi(i)=hmi(i-1)- deltahs_melt_surf(i)*ros(i)./romi*ss_fract;           % SUPERIMPOSED ICE FORMATION 
            hi(i)= hi(i-1) +deltahi_bott(i);
            
            %   EXTRA HEAT FOR ENERGY CONSERVATION
           
            if hs(i)<hs_min
                hs(i)=hs_prec_bucket(i)+(hs(i)-ZERO);          
                hs_prec_bucket(i)=0.0;
            end    
                                
%               if hs(i)<ZERO   
%                 hs(i)=hs_prec_bucket(i)+(hs(i)-ZERO);          
%                 hs_prec_bucket(i)=0.0;
%                 if hs(i)<ZERO
%                     if hmi(i)>hmi_min
%                         hmi(i)=hmi(i)+(Qs(i)*(hs(i)-ZERO)./Qmi(i)); 
%                         hs(i)=ZERO;
%                         if hmi(i)<ZERO
%                             hi(i)=hi(i)+(Qmi(i)*(hmi(i)-ZERO)./Qi_surf(i));
%                             hmi(i)=ZERO;
%                             if hi(i)<ZERO
%                                 R(i)=(Qi_surf(i)*(hi(i)-ZERO))./deltat;
%                                 hi(i)=ZERO;
%                             end
%                         end
%                     else
%                        hi(i)=hi(i)+(Qs(i)*(hs(i)-ZERO)./Qi_surf(i)); 
%                        hs(i)=ZERO;
%                        if hi(i)<ZERO
%                             R(i)=(Qi_surf(i)*(hi(i)-ZERO))./deltat;
%                             hi(i)=ZERO;
%                        end
%                     end
%                 end
%             end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SALINITY FOR MELTING SNOW
            
            if hmi(i)>hmi_min   
                Ssnowice(i)=(-deltahs_melt_surf(i)*0.20*Ssnow+hmi(i)*Ssnowice(i-1))./(hmi(i)-deltahs_melt_surf(i));
            else
                Ssnowice(i)=0.0;
            end

            if Vbr_ice(i-1)>=0.05                                                                                                   %   FLUSHING
                Q(i)=0.20*ros(i).*(-deltahs_melt_surf(i))/deltat;                                                                   %   Salinity of Q is 0      %   0.80 is refreezing, 0.20 is percolating
                Sbr_ice_star(i) = Sbr_ice(i-1) - ((Q(i)./(Vbr_ice(i-1).*ro_br(i-1)*10^3)).*(Sbr_ice(i-1) - Sbr_end))*deltat;
                Sice(i)= Sice(i-1) + (Sbr_ice_star(i)-Sbr_ice(i-1)).*Vbr_ice(i-1);                                                  % brine volume fixed in this process, it changes only as function of temperature in the next equation for next timestep                
                Vbr_ice(i)=-mu*Sice(i)./Tice(i);
                Sbr_ice(i)=Sice(i)./Vbr_ice(i);
                
                Sbr_i_star(i)=Sbr_i(i-1) - ((Q(i)./(Vbr_i(i-1).*ro_br_i(i-1)*10^3)).*(Sbr_i(i-1) - Sbr_end))*deltat;
                Si(i)= Si(i-1) + (Sbr_i_star(i)-Sbr_i(i-1)).*Vbr_i(i-1);
                Vbr_i(i)=-mu*Si(i)./(Ti(i)-273.15);
                Sbr_i(i)=Si(i)./Vbr_i(i);

                Sbr_5_star(i)=Sbr_5(i-1) - ((Q(i)./(0.05*ro_br_5(i-1)*10^3)).*(Sbr_5(i-1) - Sbr_end))*deltat;
                Sice_5(i)=Sice_5(i-1) + (Sbr_5_star(i)-Sbr_5(i-1))*0.05;               
                Tice_5(i)=-mu*Sice_5(i)/0.05;
                Sbr_5(i)=Sice_5(i)/0.05;
%                 if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))
% 
%                     hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
%                 else
%                     %hi_5(i)=ZERO;
%                     hi_5(i)=hi_5(i-1);
%                 end

                hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
                Sbr_bio_star(i)=Sbr_bio(i-1)-((Q(i)./(Vbr_bio(i-1).*ro_br_bio(i-1)*10^3)).*(Sbr_bio(i-1) - Sbr_end))*deltat;        %   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
                Sice_bio(i)=Sice_bio(i-1) + (Sbr_bio_star(i)-Sbr_bio(i-1)).*Vbr_bio(i-1);
                Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.15) + ki_ice_bio(i-1)*Tice(i))/(ki_bio_bott(i-1)+ki_ice_bott(i-1));
                Vbr_bio(i)=-mu*Sice_bio(i)./Tice_bio(i);
                Sbr_bio(i)=Sice_bio(i)./Vbr_bio(i);
                hi_bio(i)=hi_5(i);
                
                Sbr_bott_star(i)=Sbr_bott(i-1)-((Q(i)./(Vbr_bott(i-1).*ro_br_bott(i-1)*10^3)).*(Sbr_bott(i-1) - Sbr_end))*deltat;   %   DESALINATION IN THE BOTTOM LAYER
                Sice_bott(i)=Sice_bott(i-1) + (Sbr_bott_star(i)-Sbr_bott(i-1)).*Vbr_bott(i-1);
                Vbr_bott(i)=-mu*Sice_bott(i)./(Tfr-273.1499);
                Sbr_bott(i)=Sice_bott(i)./Vbr_bott(i);

            elseif Vbr_ice(i-1)<0.05                                              %   NO FLUSHING
                Q(i)=0.0;
                Sice(i)=Sice(i-1);
                Sbr_ice(i)=-Tice(i)/mu;
                Vbr_ice(i)=Sice(i)./Sbr_ice(i);
                
                Si(i)=Si(i-1);
                Sbr_i(i)=-(Ti(i)-273.15)/mu;
                Vbr_i(i)=Si(i)./Sbr_i(i);

                Sice_5(i)=Sice_5(i-1);
                Sbr_5(i)=Sice_5(i)/0.05;
                Tice_5(i)=-mu*Sbr_5(i);
%                 if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))
% 
%                     hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
%                 else
% %                     hi_5(i)=ZERO;
%                     hi_5(i)=hi_5(i-1);
%                 end
                hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
                Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.15) + ki_ice_5(i-1)*Tice_5(i))/(ki_bio_bott(i-1)+ki_ice_5(i-1));
                Sice_bio(i)=Sice_bio(i-1);                      %   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
                Sbr_bio(i)=-Tice_bio(i)/mu;
                Vbr_bio(i)=Sice_bio(i)./Sbr_bio(i);
                hi_bio(i)=hi_5(i);
                
                Sice_bott(i)=Sice_bott(i-1);
                Sbr_bott(i)=-(Tfr-273.1499)/mu;
                Vbr_bott(i)=Sice_bott(i)./Sbr_bott(i);

            end 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif ((hs(i-1)+hs_prec_bucket(i))<=hs_min)
            
            if hmi(i-1)>hmi_min
                
                %   CASE 2: MELTING SNOW ICE
                
                hs(i)=ZERO;
                T0(i)=Tfr;
                I0(i)=vis_fr*Fs(i).*expm(-kmi_av(i).*hmi(i-1));
                IM(i)=I0(i);

                if hi(i-1)>=0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i)*0.1);
                    ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
                elseif hi(i-1)<0.1
                    ISI_10(i)=IM(i).*expm(-ksi_10_av(i).*hi(i-1));
                    ISI(i)=ISI_10(i);
                end
                
                Fsurf(i)=Kmi(i).*(Tsi(i)-T0(i));
                Fbott(i)=Kice_bott(i).*(Tfr-Ti(i));
                
                delta_sublim(i)=(Fla(i)*deltat)/(romi.*Lv-Qmi(i));
                
                if Fsurf(i)>=(F(i)+I0(i))
                    deltahmi_melt_surf(i)=deltat*(F(i)+I0(i)-Fsurf(i))/Qmi(i);
                else
                    deltahmi_melt_surf(i)=0.0;
                end
                
                deltahi_melt_surf(i)=0.0;
                deltahs_melt_surf(i)=0.0;
                deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);
                
                hs(i)=ZERO+delta_bucket(i);
                hmi(i)=hmi(i-1) +deltahmi_melt_surf(i)-delta_sublim(i);
                hi(i)=  hi(i-1)+deltahi_bott(i);

                %   EXTRA HEAT FOR ENERGY CONSERVATION
                
                if hmi(i)<ZERO
                    hi(i)=hi(i)+(Qmi(i)*(hmi(i)-ZERO)./Qi_surf(i));
                    hmi(i)=ZERO;
                    if hi(i)<ZERO
                        R(i)=(Qi_surf(i)*(hi(i)-ZERO))./deltat;
                        hi(i)=ZERO;
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SALINITY FOR MELTING SNOW ICE

                Ssnowice(i)=Ssnowice(i-1);
                
                if Vbr_ice(i-1)>=0.05                                                                   %   FLUSHING
                    Q(i)=0.20*romi*(-deltahmi_melt_surf(i))/deltat;                                     %   Salinity of Q is 0      %   0.70 is refreezing, 0.20 is percolating, 0.10 is sublimating
                    Sbr_ice_star(i) = Sbr_ice(i-1) - ((Q(i)./(Vbr_ice(i-1).*ro_br(i-1)*10^3)).*(Sbr_ice(i-1) - Sbr_end))*deltat;
                    Sice(i)= Sice(i-1) + (Sbr_ice_star(i)-Sbr_ice(i-1)).*Vbr_ice(i-1);
                    Vbr_ice(i)=-mu*Sice(i)./Tice(i);
                    Sbr_ice(i)=Sice(i)./Vbr_ice(i);
                    
                    Sbr_i_star(i) = Sbr_i(i-1) - ((Q(i)./(Vbr_i(i-1).*ro_br_i(i-1)*10^3)).*(Sbr_i(i-1) - Sbr_end))*deltat;
                    Si(i)= Si(i-1) + (Sbr_i_star(i)-Sbr_i(i-1)).*Vbr_i(i-1);
                    Vbr_i(i)=-mu*Si(i)./(Ti(i)-273.15);
                    Sbr_i(i)=Si(i)./Vbr_i(i);
                    
                    Sbr_5_star(i)=Sbr_5(i-1) - ((Q(i)./(0.05*ro_br_5(i-1)*10^3)).*(Sbr_5(i-1) - Sbr_end))*deltat;
                    Sice_5(i)=Sice_5(i-1) + (Sbr_5_star(i)-Sbr_5(i-1))*0.05;
                    Tice_5(i)=-mu*Sice_5(i)/0.05;
                    Sbr_5(i)=Sice_5(i)./0.05;
%                     if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))
% 
%                         hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
%                     else
%                         % hi_5(i)=ZERO;
%                         hi_5(i)=hi_5(i-1);
%                     end

                    Sbr_bio_star(i)=Sbr_bio(i-1)-((Q(i)./(Vbr_bio(i-1).*ro_br_bio(i-1)*10^3)).*(Sbr_bio(i-1) - Sbr_end))*deltat;    %   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
                    Sice_bio(i)=Sice_bio(i-1) + (Sbr_bio_star(i)-Sbr_bio(i-1)).*Vbr_bio(i-1);
                    Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.15) + ki_ice_bio(i-1)*Tice(i))/(ki_bio_bott(i-1)+ki_ice_bott(i-1));
                    Vbr_bio(i)=-mu*Sice_bio(i)./Tice_bio(i);
                    Sbr_bio(i)=Sice_bio(i)./Vbr_bio(i);
                    hi_bio(i)=hi_5(i);
                    
                    Sbr_bott_star(i)=Sbr_bott(i-1)-((Q(i)./(Vbr_bott(i-1).*ro_br_bott(i-1)*10^3)).*(Sbr_bott(i-1) - Sbr_end))*deltat;  %   DESALINATION IN THE BOTTOM LAYER
                    Sice_bott(i)=Sice_bott(i-1) + (Sbr_bott_star(i)-Sbr_bott(i-1)).*Vbr_bott(i-1);
                    Vbr_bott(i)=-mu*Sice_bott(i)./(Tfr-273.1499);
                    Sbr_bott(i)=Sice_bott(i)./Vbr_bott(i);

                elseif Vbr_ice(i-1)<0.05                                              %   NO FLUSHING
                    Q(i)=0.0;
                    Sice(i)=Sice(i-1);
                    Sbr_ice(i)=-Tice(i)/mu;
                    Vbr_ice(i)=Sice(i)./Sbr_ice(i);
                    
                    Si(i)=Si(i-1);
                    Sbr_i(i)=-(Ti(i)-273.15)/mu;
                    Vbr_i(i)=Si(i)./Sbr_i(i);

                    Sice_5(i)=Sice_5(i-1);
                    Sbr_5(i)=Sice_5(i)/0.05;
                    Tice_5(i)=-mu*Sbr_5(i);
%                     if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))
% 
%                         hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
%                     else
%                         %hi_5(i)=ZERO;
%                         hi_5(i)=hi_5(i-1);
%                     end
                    hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
                    Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.15) + ki_ice_5(i-1)*Tice_5(i))/(ki_bio_bott(i-1)+ki_ice_5(i-1));
                    Sice_bio(i)=Sice_bio(i-1);                                      %   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
                    Sbr_bio(i)=-Tice_bio(i)/mu;
                    Vbr_bio(i)=Sice_bio(i)./Sbr_bio(i);
                    hi_bio(i)=hi_5(i);
                    
                    Sice_bott(i)=Sice_bott(i-1);
                    Sbr_bott(i)=-(Tfr-273.1499)/mu;
                    Vbr_bott(i)=Sice_bott(i)./Sbr_bott(i);


                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            elseif hmi(i-1)<=hmi_min
                
                if hi(i-1)>ZERO
                    
                    %    CASE 3: MELTING SEA ICE
                    
                    hs(i)=ZERO;
                    hmi(i)=ZERO;
                    T0(i)=Tfr;
                    if hi(i-1)>=0.1
                        ISI_10(i)=vis_fr*Fs(i).*expm(-ksi_10_av(i)*0.1);
                        ISI(i)=ISI_10(i).*expm(-ksi_av(i).*(hi(i-1)-0.1));
                    elseif hi(i-1)<0.1
                        ISI_10(i)=vis_fr*Fs(i).*expm(-ksi_10_av(i).*hi(i-1));
                        ISI(i)=ISI_10(i);
                    end

                    I0(i)=ISI_10(i);
                    IM(i)=I0(i);
                   
                    Fsurf(i)=Kice_surf(i).*(Ti(i)-T0(i));
                    Fbott(i)=Kice_bott(i).*(Tfr-Ti(i));

                    delta_sublim(i)=(Fla(i)*deltat)/(ro_sice_surf(i-1).*Lv-Qi_surf(i));
                    deltahi_bott(i)=deltat*(Fbott(i)-Fw(i)-ISI(i))./Qi_bott(i);
                    
                    if (Fsurf(i)>=(F(i)+I0(i)))
                        deltahi_melt_surf(i)=deltat*(F(i) + I0(i)-Fsurf(i))./Qi_surf(i);
                    else
                        deltahi_melt_surf(i)=0.0;
                    end
                          
                    deltahs_melt_surf(i)=0.0;
                    deltahmi_melt_surf(i)=0.0;
                    
                    hs(i)=ZERO+delta_bucket(i);
                    hmi(i)=ZERO;
                    hi(i)=  hi(i-1) + deltahi_melt_surf(i) +deltahi_bott(i) -delta_sublim(i);

                    if hi(i)<ZERO
                        R(i)=(Qi_surf(i)*(hi(i)-ZERO))./deltat;
                        hi(i)=ZERO;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % SALINITY FOR MELTING SEA ICE
                    Ssnowice(i)=0.0;
                    if hi(i)>ZERO
                        if Vbr_ice(i-1)>=0.05                                                           %   FLUSHING
                            Q(i)=0.20*ro_sice_surf(i-1)*(-deltahi_melt_surf(i))/deltat;                 %   Salinity of Q is 0      %   0.80 is refreezing, 0.20 is percolating
                            Sbr_ice_star(i) = Sbr_ice(i-1) - ((Q(i)./(Vbr_ice(i-1).*ro_br(i-1)*10^3)).*(Sbr_ice(i-1) - Sbr_end))*deltat;
                            Sice(i)= Sice(i-1) + (Sbr_ice_star(i)-Sbr_ice(i-1)).*Vbr_ice(i-1);
                            Vbr_ice(i)=-mu*Sice(i)./Tice(i);
                            Sbr_ice(i)=Sice(i)./Vbr_ice(i);
                            
                            Sbr_i_star(i) = Sbr_i(i-1) - ((Q(i)./(Vbr_i(i-1).*ro_br_i(i-1)*10^3)).*(Sbr_i(i-1) - Sbr_end))*deltat;
                            Si(i)= Si(i-1) + (Sbr_i_star(i)-Sbr_i(i-1)).*Vbr_i(i-1);
                            Vbr_i(i)=-mu*Si(i)./(Ti(i)-273.15);
                            Sbr_i(i)=Si(i)./Vbr_i(i);

                            Sbr_5_star(i)=Sbr_5(i-1) - ((Q(i)./(0.05*ro_br_5(i-1)*10^3)).*(Sbr_5(i-1) - Sbr_end))*deltat;
                            Sice_5(i)=Sice_5(i-1) + (Sbr_5_star(i)-Sbr_5(i-1))*0.05;
                            Tice_5(i)=(ki_5_bott(i-1)*(Tfr-273.15) + ki_ice_5(i-1).*Tice(i))./(ki_5_bott(i-1)+ki_ice_5(i-1));
                            Vbr_5(i)=-mu*Sice_5(i)./Tice_5(i);
                            Sbr_5(i)=Sice_5(i)./Vbr_5(i);
%                             if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))
% 
%                                 hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
%                             else
%                                 %hi_5(i)=ZERO;
%                                 hi_5(i)=hi_5(i-1);
%                             end
                            hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
                            Sbr_bio_star(i)=Sbr_bio(i-1)-((Q(i)./(Vbr_bio(i-1).*ro_br_bio(i-1)*10^3)).*(Sbr_bio(i-1) - Sbr_end))*deltat;       %   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
                            Sice_bio(i)=Sice_bio(i-1) + (Sbr_bio_star(i)-Sbr_bio(i-1)).*Vbr_bio(i-1);
                            Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.15) + ki_ice_bio(i-1)*Tice(i))/(ki_bio_bott(i-1)+ki_ice_bott(i-1));
                            Vbr_bio(i)=-mu*Sice_bio(i)./Tice_bio(i);
                            Sbr_bio(i)=Sice_bio(i)./Vbr_bio(i);
                            hi_bio(i)=hi_5(i);
                            
                            Sbr_bott_star(i)=Sbr_bott(i-1)-((Q(i)./(Vbr_bott(i-1).*ro_br_bott(i-1)*10^3)).*(Sbr_bott(i-1) - Sbr_end))*deltat;  %   DESALINATION IN THE BOTTOM LAYER
                            Sice_bott(i)=Sice_bott(i-1) + (Sbr_bott_star(i)-Sbr_bott(i-1)).*Vbr_bott(i-1);
                            Vbr_bott(i)=-mu*Sice_bott(i)./(Tfr-273.1499);
                            Sbr_bott(i)=Sice_bott(i)./Vbr_bott(i);
                            

                        elseif Vbr_ice(i-1)<0.05                                              %   NO FLUSHING
                            Q(i)=0.0;
                            Sice(i)=Sice(i-1);
                            Sbr_ice(i)=-Tice(i)/mu;
                            Vbr_ice(i)=Sice(i)./Sbr_ice(i);
                            
                            Si(i)=Si(i-1);
                            Sbr_i(i)=-(Ti(i)-273.15)/mu;
                            Vbr_i(i)=Si(i)./Sbr_i(i);
                            
                            Sice_5(i)=Sice_5(i-1);
                            Sbr_5(i)=Sice_5(i)/0.05;
                            Tice_5(i)=-mu*Sbr_5(i);
%                             if ((Tice(i)<Tice_5(i)) && (Tice_5(i)<(Tfr-273.15)))
% 
%                                 hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
%                             else
%                                 %hi_5(i)=ZERO;
%                                 hi_5(i)=hi_5(i-1);
%                             end
                            hi_5(i)=(hi(i).*ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15)))./((ki_ice_5(i-1).*(Tice(i)-Tice_5(i)))+(ki_5_bott(i-1).*(Tice_5(i)-(Tfr-273.15))));
                            Tice_bio(i)=(ki_bio_bott(i-1)*(Tfr-273.15) + ki_5_bio(i-1).*Tice_5(i))./(ki_bio_bott(i-1)+ki_5_bio(i-1));
                            Sice_bio(i)=Sice_bio(i-1);                                      %   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
                            Sbr_bio(i)=-Tice_bio(i)./mu;
                            Vbr_bio(i)=Sice_bio(i)./Sbr_bio(i);
                            hi_bio(i)=hi_5(i);
                            
                            Sice_bott(i)=Sice_bott(i-1);
                            Sbr_bott(i)=-(Tfr-273.1499)/mu;
                            Vbr_bott(i)=Sice_bott(i)./Sbr_bott(i);
                            

                        end
                    elseif hi(i)<=ZERO
                        % NO SALINITY COMPUTATION
                        Q(i)=ZERO;
                        Sice(i)=0.0;
                        Sice_bott(i)=0.0;
                        Sbr_ice(i)=0.0;
                        Vbr_ice(i)=0.0;
                        
                        Si(i)=0.0;
                        Sbr_i(i)=0.0;
                        Vbr_i(i)=0.0;

                        Sbr_5(i)=0.0;
                        Sice_5(i)=0.0;
                        Tice_5(i)=Tfr-273.15;
                        hi_5(i)=ZERO;

                        Sbr_bio(i)=0.0;
                        Sice_bio(i)=0.0;
                        Tice_bio(i)=Tfr-273.15;
                        Vbr_bio(i)=0.0;
                        hi_bio(i)=hi_5(i);
                        
                        Sice_bott(i)=0.0;
                        Sbr_bott(i)=0.0;
                        Vbr_bott(i)=0.0;
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                elseif hi(i-1)<=ZERO
                    
                    %   CASE 4: ONLY SEAWATER

                    ISI(i)=Fs(i).*((infra_fr*expm(-k_ocean_red*h_mix)) + ((1-infra_fr)*expm(-k_ocean_vis*h_mix)));        % Penetrating solar radiation in sea water (type II, Jerlov, 1968)
                    I0(i)=ISI(i);
                    
                    Fsurf(i)=0.0;
                    Fbott(i)=0.0;

                    IM(i)=ISI(i);
                    ISI_10(i)=ISI(i);

                    deltahs_melt_surf(i)=0.0;
                    deltahmi_melt_surf(i)=0.0;
                    deltahi_melt_surf(i)=0.0;
                    deltahi_bott(i)=0.0;
                    delta_bucket(i)=0.0;
                    hs_prec_bucket(i)=0.0;

                    hs(i)=ZERO;
                    hi(i)=ZERO;
                    hmi(i)=ZERO;

                    Ri(i)=0.0;

                    % NO SALINITY COMPUTATION

                    Q(i)=ZERO;
                    Sice(i)=0.0;
                    Sice_bott(i)=0.0;
                    Sbr_ice(i)=0.0;
                    Vbr_ice(i)=0.0;
                    Ssnowice(i)=0.0;
                    
                    Si(i)=0.0;
                    Sbr_i(i)=0.0;
                    Vbr_i(i)=0.0;

                    Sbr_5(i)=0.0;
                    Sice_5(i)=0.0;
                    Tice_5(i)=Tfr-273.15;
                    hi_5(i)=ZERO;

                    Sbr_bio(i)=0;
                    Sice_bio(i)=0.0;
                    Tice_bio(i)=Tfr-273.15;
                    Vbr_bio(i)=0.0;
                    hi_bio(i)=hi_5(i);
                    
                    Sice_bott(i)=0.0;
                    Sbr_bott(i)=0.0;
                    Vbr_bott(i)=0.0;

                end
            end
        end
    else
        'ERROR MELTING SEASON'
    end
    
    
    % PROPERTIES DEPENDING ON T,S

    Tice_av(i)=Ti(i)-273.15;                        %   =====================>>>>>>>>>>>>>>>> TO BFM
    Sice_av(i)=(Sice_bott(i)+Sice(i))/2;     %   =====================>>>>>>>>>>>>>>>> TO BFM

    if Tice_av(i) > -2.000
        F1(i)= -4.1221*10^(-2) - 1.8407*10*Tice_av(i) + 5.4802*10^(-1)*Tice_av(i).^2 + 2.1454*10^(-1)*Tice_av(i).^3;
        F2(i)= 9.0312*10^(-2) -1.6111*10^(-2)*Tice_av(i) +1.2291*10^(-4)*Tice_av(i).^2 +1.3603*10^(-4)*Tice_av(i).^3;
    else
        F1(i)= -4.723 -2.245*10*Tice_av(i) - 6.397*10^(-1)*Tice_av(i).^2 -1.074*10^(-2)*Tice_av(i).^3;
        F2(i)= 8.903*10^(-2) -1.763*10^(-2)*Tice_av(i) -5.330*10^(-4)*Tice_av(i).^2 -8.801*10^(-4)*Tice_av(i).^3;
    end

    roi_surf(i)=0.917-1.403*10^(-4)*Tice(i);                                            % density of pure ice at the surface
    ro_sice_surf(i)=0.905;%(1-Va)*((roi_surf(i).*F1(i))./(F1(i)-roi_surf(i).*Sice(i).*F2(i))); % bulk density of sea ice at the surface as function of T,S;
    ro_sice_bott(i)=0.885;%(1-Va)*((roi_bott*F1(i))./(F1(i)-roi_bott*Sice_bott(i).*F2(i)));    % bulk density of sea ice at the bottm (density of pure ice at the bottom is fixed at 0.900)
    roi_av(i)=0.917-1.403*10^(-4)*Tice_av(i);                                           % density of pure ice of the layer as function of Tice_av in g/cm^3
    ro_sice_bulk(i)=(1-Va)*((roi_av(i).*F1(i))./(F1(i)-roi_av(i).*Sice_av(i).*F2(i)));  % bulk density of ice as function of Tice,Sice in g/cm^3

    ro_br(i)=1+8*10^(-4)*Sbr_ice(i);                                                    % brines density as function of brines salinity in g/cm^3
    ro_br_i(i)=1+8*10^(-4)*Sbr_i(i);   
    ro_br_5(i)=1+8*10^(-4)*Sbr_5(i);
    ro_br_bio(i)=1+8*10^(-4)*Sbr_bio(i);
    ro_br_bott(i)=1+8*10^(-4)*Sbr_bott(i);

    k0i(i) = 418.6*(5.35*10^(-3) - 2.568*10^(-5)*Tice_av(i));                           % thermal conductivity of pure ice as function of T
    kb(i) = 418.6*(1.25*10^(-3) + 3.0*10^(-5)*Tice_av(i));                              % thermal conductivity of pure brines as function of T

    k0i_5_bio(i)=418.6*(5.35*10^(-3) - 2.568*10^(-5)*((Tice_5(i)+Tice_bio(i))/2));
    kb_5_bio(i)=418.6*(1.25*10^(-3) + 3.0*10^(-5)*((Tice_5(i)+Tice_bio(i))/2));
    ki_5_bio(i)= (1-Va-0.05)*k0i_5_bio(i) + 0.005*kb_5_bio(i);

    k0i_ice_5(i)=418.6*(5.35*10^(-3) - 2.568*10^(-5)*((Tice(i)+Tice_5(i))/2));
    kb_ice_5(i)=418.6*(1.25*10^(-3) + 3.0*10^(-5)*((Tice(i)+Tice_5(i))/2));
    ki_ice_5(i)=(1-Va-Vbr_ice(i)).*k0i_ice_5(i) + Vbr_ice(i).*kb_ice_5(i);

    k0i_ice_bott(i)=418.6*(5.35*10^(-3) - 2.568*10^(-5)*((Tice(i)+(Tfr-273.15))/2));
    kb_ice_bott(i)=418.6*(1.25*10^(-3) + 3.0*10^(-5)*((Tice(i)+(Tfr-273.15))/2));
    ki_ice_bott(i)=(1-Va-Vbr_ice(i)).*k0i_ice_bott(i) + Vbr_ice(i).*kb_ice_bott(i);

    k0i_5_bott(i)=418.6*(5.35*10^(-3) - 2.568*10^(-5)*((Tice_5(i)+(Tfr-273.15))/2));
    kb_5_bott(i)=418.6*(1.25*10^(-3) + 3.0*10^(-5)*((Tice_5(i)+(Tfr-273.15))/2));
    ki_5_bott(i)=(1-Va-0.05)*k0i_5_bott(i) + 0.05*kb_5_bott(i);

    k0i_ice_bio(i)=418.6*(5.35*10^(-3) - 2.568*10^(-5)*((Tice(i)+Tice_bio(i))/2));
    kb_ice_bio(i)=418.6*(1.25*10^(-3) + 3.0*10^(-5)*((Tice(i)+Tice_bio(i))/2));
    ki_ice_bio(i)= (1-Va-Vbr_ice(i)).*k0i_ice_bio(i) + Vbr_ice(i).*kb_ice_bio(i);

    k0i_bio_bott(i)=418.6*(5.35*10^(-3) - 2.568*10^(-5)*((Tice_bio(i)+(Tfr-273.15))/2));
    kb_bio_bott(i)=418.6*(1.25*10^(-3) + 3.0*10^(-5)*((Tice_bio(i)+(Tfr-273.15))/2));
    ki_bio_bott(i)= (1-Va-Vbr_bio(i)).*k0i_bio_bott(i) + Vbr_bio(i).*kb_bio_bott(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   COMPUTATION OF MEAN PROPERTIES OF THE BIOLOGICAL ACTIVE PART OF THE SEA ICE SISTEM FOR BFM RUN

    if hi(i)>hi_min
        if hi(i)>=0.1
            ISI_layer(i)=I0(i).*expm(-ksi_av(i).*((hi(i)-0.1)/2));
            ISI_bio(i)=I0(i).*expm(-ksi_av(i).*(hi(i)-0.1- hi_5(i)/2));

        elseif hi(i)<0.1
            ISI_layer(i)=IM(i).*expm(-ksi_10_av(i).*(hi(i)/2));
            ISI_bio(i)=IM(i).*expm(-ksi_10_av(i).*(hi(i)-hi_5(i)/2));

        end

    else
        ISI_layer(i)=0.0;
        ISI_bio(i)=0.0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END RUNNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENT OUT WHEN TESTING MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BFM INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BFM = [ESIM]

% EHB = [hi_bio]				sea ice biotic layer thickness (BAL) 
% EHI = [hi+hmi]				total ice thickness 
% EVB = [Vbr_bio]				simulated brine volume between 0.05 and the bottom brine volume 
% ETB = [Tice_bio]				simulated Temperature of ice/brines between Ti_5 and the bottom freezing temperature 
% ESB = [Sbr_bio]				simulated brines salinity between Sbr_5 and the bottom brine salinity 
% EIB = [ISI_bio]				simulated irradiance at the middle point in the biological active system 
% ESI = [Sice_bio]				simulated bulk salinity at the middle point in the BAL
% EDH = [deltahi_bott(i)./deltat]		simulated ice growth/melt velocity 
% EDS = [-deltahs_melt_surf(i)./deltat]		simulated snow melt velocity (from negative in ESIM to positive when enter BFM) 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!! BE AWARE: TO RUN AFTER ESIM MODEL RUNS!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  for i=1:nmax
%     if hi_bio(i)>ZERO
%         EHB_i(i)=hi_bio(i);                          	% simulated ice thickness of the biological active layer
%         EHI_i(i)=hi(i)+hmi(i);                        	% total ice thickness
% %         if EHI_i(i)< EHB_i(i)
% %             EHI_i(i)=EHB_i(i);
% %         end
%         if EHB_i(i)> EHI_i(i)
%             EHB_i(i)=EHI_i(i);
%         end
%         EVB_i(i)=Vbr_bio(i);                         	% simulated brine volume between 0.05 and the bottom brine volume
%         ETB_i(i)=Tice_bio(i);                        	% simulated Temperature of ice/brines between Ti_5 and the bottom freezing temperature
%         ESB_i(i)=Sbr_bio(i);                         	% simulated brines salinity between Sbr_5 and the bottom brine salinity
%         if hi_bio(i)>=0.1
%             ISI_bio(i)=I0(i).*expm(-ksi_av(i).*(hi(i)-0.1- hi_bio(i)/2));
%         elseif hi(i)<0.1
%             ISI_bio(i)=IM(i).*expm(-ksi_10_av(i).*(hi(i)-hi_bio(i)/2));
%         end
%         EIB_i(i)=ISI_bio(i);                       	% simulated irradiance at the middle point in the biological active system
%         ESI_i(i)=Sice_bio(i);                        	% simulated bulk salinity at the middle point in the biological active system
%         EDH_i(i)=deltahi_bott(i)./deltat;            	% simulated ice growth/melt velocity
%         EDS_i(i)=-deltahs_melt_surf(i)./deltat;      	% simulated snow melt velocity (from negative in ESIM to positive when enter BFM)
%         EICE_i(i)=1;
%         EHI=hi+hmi;                                	% total ice thickness
%     else
%         EHB_i(i)=0.0;
%         EHI_i(i)=hi(i)+hmi(i);
%         EVB_i(i)=0.0;                               
%         ETB_i(i)=Tfr-273.15;                               
%         ESB_i(i)=0.0;                            
%         EIB_i(i)=Fs(i); 
%         ESI_i(i)=Sw(i);
%         EDH_i(i)=0.0;
%         EDS_i(i)=0.0;
%         EICE_i(i)=0;
%     end   
%  end
%  for i=1:nmax 
%      T_av_i(i)=Ti(i)-273.15;
%      Sbr_av_i(i)=Sbr_i(i);
%      Sbk_av_i(i)=Si(i);
%  end
%  
%  
%  EHB=EHB_i;
%  EHI=EHI_i;
%  EVB=EVB_i;                               
%  ETB=ETB_i;                               
%  ESB=ESB_i;                            
%  EIB=EIB_i; 
%  ESI=ESI_i;
%  EDH=EDH_i;
%  EDS=EDS_i;
%  EICE=EICE_i;
%  ETW=Tmix-273.15;
%  ESW=Sw;
%  EWIND=Ua;
%  EIW=Fs;                           % Irradiance at the water surface   
%  ESH=hs+hs_prec_bucket;            % snow thickness
%  EMI=hmi;                          % intermediate layer thickness 
% 
%  T_av=T_av_i;
%  Sbr_av=Sbr_av_i;
%  Sbk_av=Sbk_av_i;
%  
% 
% ESH(5476:5837)=[];			%These are specific for Antarctic test case study site. To be adjusted by the user as needed.
% EHI(5476:5837)=[];
% EMI(5476:5837)=[];
% 
% EICE(5476:5837)=[];
% EHB(5476:5837)=[];
% EVB(5476:5837)=[];
% ETB(5476:5837)=[];
% ESB(5476:5837)=[];
% EIB(5476:5837)=[];
% ESI(5476:5837)=[];
% EDH(5476:5837)=[];
% EDS(5476:5837)=[];
% 
% ETW(5476:5837)=[];
% ESW(5476:5837)=[];
% EWIND(5476:5837)=[];
% EIW(5476:5837)=[];
% 
% T_av(5476:5837)=[];
% Sbr_av(5476:5837)=[];
% Sbk_av(5476:5837)=[];
% 
% ESH=reshape(ESH,15,365);
% EHI=reshape(EHI,15,365);
% EMI=reshape(EMI,15,365);
% 
% EICE=reshape(EICE,15,365);
% EHB=reshape(EHB,15,365);
% EVB=reshape(EVB,15,365);
% ETB=reshape(ETB,15,365);
% ESB=reshape(ESB,15,365);
% EIB=reshape(EIB,15,365);
% ESI=reshape(ESI,15,365);
% EDH=reshape(EDH,15,365);
% EDS=reshape(EDS,15,365);
% 
% ETW=reshape(ETW,15,365);
% ESW=reshape(ESW,15,365);
% EWIND=reshape(EWIND,15,365);
% EIW=reshape(EIW,15,365);
% 
% T_av=reshape(T_av,15,365);
% Sbr_av=reshape(Sbr_av,15,365);
% Sbk_av=reshape(Sbk_av,15,365);
% 
% 
% ESH=mean(ESH);
% EHI=mean(EHI);
% EMI=mean(EMI);
% 
% EICE=mean(EICE);
% EHB=mean(EHB);
% EVB=mean(EVB);
% ETB=mean(ETB);
% ESB=mean(ESB);
% EIB=mean(EIB);
% ESI=mean(ESI);
% EDH=mean(EDH);
% EDS=mean(EDS);
% 
% ETW=mean(ETW);
% ESW=mean(ESW);
% EWIND=mean(EWIND);
% EIW=mean(EIW);
% 
% T_av=mean(T_av);
% Sbr_av=mean(Sbr_av);
% Sbk_av=mean(Sbk_av);
% 
% ESH=ESH';
% EHI=EHI';
% EMI=EMI';
% 
% EICE=EICE';
% EIB=EIB';
% EVB=EVB';
% ETB=ETB';
% ESB=ESB';
% ESI=ESI';
% EHB=EHB';
% EDH=EDH';
% EDS=EDS';
% 
% ETW=ETW';
% ESW=ESW';
% EWIND=EWIND';
% EIW=EIW';
% 
% T_av=T_av';
% Sbr_av=Sbr_av';
% Sbk_av=Sbk_av';
% 
% 
% day=736696:1:737060; % 01-01-17 : 31:12:17 (365d)		This is specific for ANTARCTIC test case study site. To be adjusted by the user as needed.
% day=day';
% day=datestr(day,31);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BFM_input=[day EICE EHB EVB ETB ESB EIB ESI];
% ESIM_forcing_num=[EICE EHB EVB ETB ESB EIB ESI EDH EDS];
% PEL_forcing_num=[EICE ETW ESW EWIND EIW];
% 
% ESIM_forcing_str=num2str(ESIM_forcing_num);
% Pel_forcing_str=num2str(PEL_forcing_num);
% 
% time=num2str(day);
% 
% BFM_input_seaice_ANT=[time ESIM_forcing_str];
% BFM_input_pelagic_ANT=[time Pel_forcing_str];
% 
% dlmwrite('BFM_input_seaice_ANT.dat', BFM_input_seaice_ANT, 'delimiter', '%s',  'precision', 6)
% dlmwrite('BFM_input_pelagic_ANT.dat', BFM_input_pelagic_ANT, 'delimiter', '%s',  'precision', 6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%it 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIST OF REFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ABOUT SIM
%   Tedesco, L., Vichi, M., Haapala, J., Stipa, T., 2008. An enhanced sea
%   ice thermodynamic model applied to the Baltic Sea, Boreal Environmental Research, in
%   press

%   ABOUT PARAMETERIZATIONS
%   Assur, A., 1958. Composition of sea ice and its tensile strenght. In: Arctic Sea Ice , Publ. Nat. Acad. Sci.
%   U.S.A., 598: 106-138.
%   Bitz and Lipscomb, 1999
%   Cox and Weeks,1988
%   Guest, P.S., 1997. Surface longwave radiation conditions in the Eastern Weddel Sea during winter, Journal of
%   Geophysical Research, 103, C13: 30761--30772.
%   Lepparanta, M, 1983. A Growth Model for Black Ice, Snow Ice and Snow Thickness in Subarctic Basins. Nordic Hydrology,
%   59--70.
%   Maykut, G.A., 1982. Large-scale heat exchange and ice-production in the central Arctic. Journal of Geophysical
%   Research, 87, C10: 7971--7984.
%   Omstedt, A. and Wettlaufer, J.S., 1992. Ice growth and oceanic heat flux: Modelingnts. Journal of
%   Geophysical Research 97, C6. Doi: 10.1029/92JC00815, issn: 0148-0227.
%   Untersteiner, 1961
%   Vancoppenolle, M., Bitz, C.M., Fichefet, T., 2007. Summer landfast sea ice desalination at Point Barrow,
%   Alaska: Modeling and observations. Journal of Geophysical Research, 112, C04022. Doi: 10.1029/2006JC003493.
%
%   ABOUT OBSERVATIONS DATA
%   Sein??, A., Peltola, J. , 1991. J????talven kestoaika ja kiintoj????n paksuustilastoja merialueilla 1961-1990 /
%   Duration of the ice season and statistics of fast ice thickness along the Finnish coast 1961-1990.
%   Finnish Marine Research, No. 258.



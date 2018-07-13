%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESIM3: STANDARD RUN FOR BFM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riesna R Audh

% Revised script from Letizia Tedesco (2007)

% Standard run of the ESIM2 for the Antarctic to produce the BFM-SI input file

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL PARAMETERS AND CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !! SET BY USER !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M O D E L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



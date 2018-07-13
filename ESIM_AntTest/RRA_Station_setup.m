%% ANTARCTIC TEST CASE ALONG 30S %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  !!! TO BE SET UP BY THE USER !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=2700;                               % = deltat = time step (s)
nmax=5837;                             % maximum number of time steps (-)
hs_prec_min=0.02;                      % minimum snow precipitation (m) to activate the model                   !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hs_min=0.02;                           % minimum snow thickness (m) to activate the model                       !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hmi_min=0.05;                          % minimum snow ice/supeimposed ice thickness (m) to activate the model   !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
hi_min=0.05;                           % minimum sea ice thickness (m)  to activate the model                   !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
ZERO=10^-5;                            % ZERO OF THE MODEL                                                      !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
Fwater=9.52;                            % = Fw = oceanic heat fluxes (W m-2)                                     !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS PARAMETER!!!
Tfreez=271.35;                       % = Tfr = -1.8 C seawater freezing temperature (K)                       !!!BE AWARE: THIS STRONGLY DEPENDS ON THE LOCATION!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

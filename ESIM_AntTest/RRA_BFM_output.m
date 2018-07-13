%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT FOR BFM INPUT based on Kobbefjord test case simulation                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!! BE AWARE: TO RUN AFTER ESIM MODEL RUNS!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:nmax
    if hi_bio(i)>ZERO
        EHB_i(i)=hi_bio(i);                          	% simulated ice thickness of the biological active layer
        EHI_i(i)=hi(i)+hmi(i);                        	% total ice thickness
%         if EHI_i(i)< EHB_i(i)
%             EHI_i(i)=EHB_i(i);
%         end
        if EHB_i(i)> EHI_i(i)
            EHB_i(i)=EHI_i(i);
        end
        EVB_i(i)=Vbr_bio(i);                         	% simulated brine volume between 0.05 and the bottom brine volume
        ETB_i(i)=Tice_bio(i);                        	% simulated Temperature of ice/brines between Ti_5 and the bottom freezing temperature
        ESB_i(i)=Sbr_bio(i);                         	% simulated brines salinity between Sbr_5 and the bottom brine salinity
        if hi_bio(i)>=0.1
            ISI_bio(i)=I0(i).*expm(-ksi_av(i).*(hi(i)-0.1- hi_bio(i)/2));
        elseif hi(i)<0.1
            ISI_bio(i)=IM(i).*expm(-ksi_10_av(i).*(hi(i)-hi_bio(i)/2));
        end
        EIB_i(i)=ISI_bio(i);                       	% simulated irradiance at the middle point in the biological active system
        ESI_i(i)=Sice_bio(i);                        	% simulated bulk salinity at the middle point in the biological active system
        EDH_i(i)=deltahi_bott(i)./deltat;            	% simulated ice growth/melt velocity
        EDS_i(i)=-deltahs_melt_surf(i)./deltat;      	% simulated snow melt velocity (from negative in ESIM to positive when enter BFM)
        EICE_i(i)=1;
        EHI=hi+hmi;                                	% total ice thickness
    else
        EHB_i(i)=0.0;
        EHI_i(i)=hi(i)+hmi(i);
        EVB_i(i)=0.0;                               
        ETB_i(i)=Tfr-273.15;                               
        ESB_i(i)=0.0;                            
        EIB_i(i)=Fs(i); 
        ESI_i(i)=Sw(i);
        EDH_i(i)=0.0;
        EDS_i(i)=0.0;
        EICE_i(i)=0;
    end   
 end
 for i=1:nmax 
     T_av_i(i)=Ti(i)-273.15;
     Sbr_av_i(i)=Sbr_i(i);
     Sbk_av_i(i)=Si(i);
 end
 
 
 EHB=EHB_i;
 EHI=EHI_i;
 EVB=EVB_i;                               
 ETB=ETB_i;                               
 ESB=ESB_i;                            
 EIB=EIB_i; 
 ESI=ESI_i;
 EDH=EDH_i;
 EDS=EDS_i;
 EICE=EICE_i;
 ETW=Tmix-273.15;
 ESW=Sw;
 EWIND=Ua;
 EIW=Fs;                           % Irradiance at the water surface   
 ESH=hs+hs_prec_bucket;            % snow thickness
 EMI=hmi;                          % intermediate layer thickness 

 T_av=T_av_i;
 Sbr_av=Sbr_av_i;
 Sbk_av=Sbk_av_i;
 

ESH(5476:5837)=[];			%These are specific for Antarctic test case study site. To be adjusted by the user as needed.
EHI(5476:5837)=[];
EMI(5476:5837)=[];

EICE(5476:5837)=[];
EHB(5476:5837)=[];
EVB(5476:5837)=[];
ETB(5476:5837)=[];
ESB(5476:5837)=[];
EIB(5476:5837)=[];
ESI(5476:5837)=[];
EDH(5476:5837)=[];
EDS(5476:5837)=[];

ETW(5476:5837)=[];
ESW(5476:5837)=[];
EWIND(5476:5837)=[];
EIW(5476:5837)=[];

T_av(5476:5837)=[];
Sbr_av(5476:5837)=[];
Sbk_av(5476:5837)=[];

ESH=reshape(ESH,15,365);
EHI=reshape(EHI,15,365);
EMI=reshape(EMI,15,365);

EICE=reshape(EICE,15,365);
EHB=reshape(EHB,15,365);
EVB=reshape(EVB,15,365);
ETB=reshape(ETB,15,365);
ESB=reshape(ESB,15,365);
EIB=reshape(EIB,15,365);
ESI=reshape(ESI,15,365);
EDH=reshape(EDH,15,365);
EDS=reshape(EDS,15,365);

ETW=reshape(ETW,15,365);
ESW=reshape(ESW,15,365);
EWIND=reshape(EWIND,15,365);
EIW=reshape(EIW,15,365);

T_av=reshape(T_av,15,365);
Sbr_av=reshape(Sbr_av,15,365);
Sbk_av=reshape(Sbk_av,15,365);


ESH=mean(ESH);
EHI=mean(EHI);
EMI=mean(EMI);

EICE=mean(EICE);
EHB=mean(EHB);
EVB=mean(EVB);
ETB=mean(ETB);
ESB=mean(ESB);
EIB=mean(EIB);
ESI=mean(ESI);
EDH=mean(EDH);
EDS=mean(EDS);

ETW=mean(ETW);
ESW=mean(ESW);
EWIND=mean(EWIND);
EIW=mean(EIW);

T_av=mean(T_av);
Sbr_av=mean(Sbr_av);
Sbk_av=mean(Sbk_av);

ESH=ESH';
EHI=EHI';
EMI=EMI';

EICE=EICE';
EIB=EIB';
EVB=EVB';
ETB=ETB';
ESB=ESB';
ESI=ESI';
EHB=EHB';
EDH=EDH';
EDS=EDS';

ETW=ETW';
ESW=ESW';
EWIND=EWIND';
EIW=EIW';

T_av=T_av';
Sbr_av=Sbr_av';
Sbk_av=Sbk_av';


day=736696:1:737060; % 01-01-17 : 31:12:17 (365d)		This is specific for Kobbefjord test case study site. To be adjusted by the user as needed.
day=day';
day=datestr(day,31);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BFM_input=[day EICE EHB EVB ETB ESB EIB ESI];
ESIM_forcing_num=[EICE EHB EVB ETB ESB EIB ESI EDH EDS];
PEL_forcing_num=[EICE ETW ESW EWIND EIW];

ESIM_forcing_str=num2str(ESIM_forcing_num);
Pel_forcing_str=num2str(PEL_forcing_num);

time=num2str(day);

BFM_input_seaice_ANT=[time ESIM_forcing_str];
BFM_input_pelagic_ANT=[time Pel_forcing_str];

dlmwrite('BFM_input_seaice_ANT.dat', BFM_input_seaice_ANT, 'delimiter', '%s',  'precision', 6)
dlmwrite('BFM_input_pelagic_ANT.dat', BFM_input_pelagic_ANT, 'delimiter', '%s',  'precision', 6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

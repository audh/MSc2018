%% ESIM OUTPUT DISPLAY 365 POINTS %%
% ESIM OUTPUT VARIABLES TRANSFORMED FROM 5837 TO 365 POINTS TO REDUUCE
% NOISE IN PLOTS
% BFM_OUTPUT.M SCRIPT USED TO CREATE THIS SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT FOR BFM INPUT based on ANTARCTIC test case simulation                   
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
EHB(5476:5837)=[];
EHB=reshape(EHB,15,365);
EHB=mean(EHB);
EHB=EHB';
dlmwrite('NCEP85_BAL',EHB)

HS=hs+hs_prec_bucket;
HS=HS';
HS(5476:5837)=[];
HS=reshape(HS,15,365);
HS=mean(HS);
HS=HS';

HMI=hmi;
HMI=HMI';
HMI(5476:5837)=[];
HMI=reshape(HMI,15,365);
HMI=mean(HMI);
HMI=HMI';

SI=hmi-hi;
SI=SI';
SI(5476:5837)=[];
SI=reshape(SI,15,365);
SI=mean(SI);
SI=SI';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

plot(HS,'Color',[.31 .25 .74],'LineWidth',1)
hold on

plot(-HMI,'c','LineWidth',1) 
hold on
plot(-EHB,'k','Linewidth',1)    %-abs(hi_bio) is plotting only the BAL values below the surface
hold on 

plot(SI,'b','LineWidth',1)
hold on

title('Thickness and cover 1|ERA 5| 15.5 W/m^2','FontSize',20,'FontWeight','bold')

%title('Thickness and cover 1|5400|8.5|365')
legend('snow cover','meteoric ice','BAL','sea ice','FontSize',12,'FontWeight','bold',3)

xlim([150 365])
set(gca,'XTick',150:30.4167:365)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',20,'FontWeight','bold')
% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-0.9 0.35])
ylabel('Thickness (m)','FontSize',22,'FontWeight','bold')
set(gca,'YTickLabel',{'-0.8','-0.6','-0.4','-0.2','0.0','0.2','0.4'},'FontSize',20,'FontWeight','bold')
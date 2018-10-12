%% average out to 365 points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT FOR BFM INPUT based on ANTARCTIC test case simulation                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NCEP: RUN AFTER NCEP INTERPOLATION AND ESIM3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SI=hmi-hi;
% SI=SI';
% SI(5476:5837)=[];
% SI=reshape(SI,15,365);
% SI=mean(SI);
% SI=SI';
% SIT=['SI_NCEP365.txt'];
% dlmwrite(SIT,SI);
% 
% Ta=Ta;
% Ta(5476:5837)=[];
% Ta=reshape(Ta,15,365);
% Ta=mean(Ta);
% Ta=Ta';
% Ta2=['Ta_NCEP365.txt'];
% dlmwrite(Ta2,Ta);
% 
% dswrf=Fsd_cloud;
% dswrf(5476:5837)=[];
% dswrf=reshape(dswrf,15,365);
% dswrf=mean(dswrf);
% dswrf=dswrf';
% swr=['dswrf_NCEP365.txt'];
% dlmwrite(swr,dswrf);
% 
% PPT=P_rate;
% PPT(5476:5837)=[];
% PPT=reshape(PPT,15,365);
% PPT=mean(PPT);
% PPT=PPT';
% TP=['PPT_NCEP365.txt'];
% dlmwrite(TP,PPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERA 5: RUN AFTER ERA5 INTERPOLATION AND ESIM3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SI=hmi-hi;
% SI=SI';
% SI(5476:5837)=[];
% SI=reshape(SI,15,365);
% SI=mean(SI);
% SI=SI';
% SIT=['SI_ERA5365.txt'];
% dlmwrite(SIT,SI);
% 
% Ta=Ta;
% Ta(5476:5837)=[];
% Ta=reshape(Ta,15,365);
% Ta=mean(Ta);
% Ta=Ta';
% Ta2=['Ta_ERA5365.txt'];
% dlmwrite(Ta2,Ta);
% 
% dswrf=Fsd_cloud;
% dswrf(5476:5837)=[];
% dswrf=reshape(dswrf,15,365);
% dswrf=mean(dswrf);
% dswrf=dswrf';
% swr=['dswrf_ERA5365.txt'];
% dlmwrite(swr,dswrf);
% 
% PPT=P_rate;
% PPT(5476:5837)=[];
% PPT=reshape(PPT,15,365);
% PPT=mean(PPT);
% PPT=PPT';
% TP=['PPT_ERA5365.txt'];
% dlmwrite(TP,PPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Ta_ERA5365.txt
load Ta_NCEP365.txt
% load EI_BAL.txt
% load NCEP.txt
% 

plot(Ta_ERA5365,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.5) %ncep
hold on

plot(Ta_NCEP365,'Color',[0.30, 0.70, 0.10],'LineWidth',1.5) %ncep
hold on

% plot(NCEP,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.5) %ncep
% hold on

% plot(EI_BAL,'Color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5)
% hold on
% 
% plot(E5_BAL,'Color',[0, 0.75, 0.75],'LineWidth',1.5) %E5
hold on

title('S-Index | Air Temperature','FontSize',25,'FontWeight','bold')

% title('BAL Thickness comparison|8.5 W/m^2','FontSize',25,'FontWeight','bold')
% legend('NCEP','ERA INTERIM','ERA 5', 'FontSize',20,'FontWeight','bold',2)

xlim([0 365])
set(gca,'XTick',0:30:365)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',23,'FontWeight','bold')

% set(gca,'YTick',-0.6:0.2:0.4)
ylim([0 0.7])
ylabel('thickness[m]','FontSize',25,'FontWeight','bold')
set(gca,'YTickLabel',{'0.0','0.05','0.10','0.15','0.20'},'FontSize',23,'FontWeight','bold')

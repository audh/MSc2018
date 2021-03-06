%% ice concentration extraction and plotting

% file='/home/riesna/Desktop/working_nc/icec.sfc.gauss.2017_daily.nc' ;
% ncdisp(file)
% 
% lon=ncread(file,'lon');
% lat=ncread(file,'lat');
% time=ncread(file,'time');
% 
% b2= 'NCEP_SICdaily';
% b3= 'ERA5_SIC';
% b4= 'EI_SIC';
% 
% SIC= ncread(file,'icec',[17 81 1],[1 1 inf],[1 1 1]);
% NCEP_SIC= reshape(squeeze(SIC(:,1,:,:)),365,1);
% f2=[b2,'.txt'];
% dlmwrite(f2,NCEP_SIC); 

% SIC= ncread(file,'ci',[121 616 1],[1 1 inf],[1 1 1]);
% ERA5= reshape(squeeze(SIC(:,1,:,:)),1460,1);
% f2=[b3,'.txt'];
% dlmwrite(f2,ERA5); 

% SIC= ncread(file,'ci',[21 36 1],[1 1 inf],[1 1 1]);
% EI= reshape(squeeze(SIC(:,1,:,:)),120,1);
% f2=[b4,'_Nov','.txt'];
% dlmwrite(f2,EI); 

% BAL=hi_bio;
% BAL=BAL';
% BAL(5476:5837)=[];
% BAL=reshape(BAL,15,365);
% BAL=mean(BAL);
% BAL=BAL';
% o=['ERA5125_BAL.txt'];
% dlmwrite(o,BAL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load E5.txt
% load EI.txt
% load NCEP.txt
% 
% figure
% 
% plot(E5,'Color',[0, 0.75, 0.75],'LineWidth',2)
% hold on
% % 
% % plot(EI,'Color',[0.4940, 0.1840, 0.5560],'LineWidth',2) 
% % hold on
% % 
% % plot(NCEP,'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)    
% % hold on 
% 
% title('ERA 5 Sea Ice Concentration','FontSize',20,'FontWeight','bold')
% % legend('ERA 5','Era Interim','NCEP','FontSize',10,'FontWeight','bold',3)
% 
% xlim([486 1460])
% set(gca,'XTick',486:121.6667:1400)
% % set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',25,'FontWeight','bold')
% set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',20,'FontWeight','bold')
% % set(gca,'YTick',-0.6:0.2:0.4)
% % ylim([0 102])
% ylabel('Sea Ice Concentration[%]','FontSize',25,'FontWeight','bold')
% % set(gca,'YTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2'},'FontSize',10,'FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ERA585_BAL.txt
load NCEP85_BAL.txt
% load ASMR_SIC.txt
% load SI35_ERA5365.txt

% NCEP2_SIC=NCEP_SIC;
% NCEP2_SIC=NCEP2_SIC';
% % NCEP2_SIC(5476:5837)=[];
% NCEP2_SIC=reshape(NCEP2_SIC,4,365);
% NCEP2_SIC=mean(NCEP2_SIC);
% NCEP2_SIC=NCEP2_SIC';
% v=['NCEP2_SIC365.txt'];
% dlmwrite(v,NCEP2_SIC);
% 
% BAL=E5_SIC;
% BAL=BAL';
% %E52_SIC(5476:5837)=[];
% BAL=reshape(BAL,4,365);
% BAL=mean(BAL);
% BAL=BAL';
% n=['E5_SIC365.txt'];
% dlmwrite(n,BAL);

% load NCEP.txt
% 
plot(-abs(NCEP85_BAL),'Color',[0.900, 0.50, 0.480],'LineWidth',2) %ncep
hold on

% plot(ASMR_SIC,'Color',[0.35, 0.40, 0.690],'LineWidth',2)
% hold on
% 
plot(-abs(ERA585_BAL),'Color',[0, 0.75, 0.75],'LineWidth',2) %E5
hold on

title('Sea Ice Extent at 30°E,63°S','FontSize',25,'FontWeight','bold')
legend('NCEP','ASMR','ERA5', 'FontSize',20,'FontWeight','bold',3)

xlim([150 365])
% set(gca,'XTick',150:10:195)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% set(gca,'XTickLabel',{'jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',23,'FontWeight','bold')
% set(gca,'XTickLabel',{'01-06-17','10-06-17','20-06-17','01-07-17','10-07-17'},'FontSize',23,'FontWeight','bold')

% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-1 1])
ylabel('Concentration[%]','FontSize',25,'FontWeight','bold')
% set(gca,'YTickLabel',{'0','10','20','30','40','50','60','70','80','90','100'},'FontSize',23,'FontWeight','bold')

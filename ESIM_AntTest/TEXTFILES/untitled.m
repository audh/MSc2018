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

load SI125_NCEP365.txt
load SI125_ERA5365.txt
load ASMR_SIC.txt
% load SI35_ERA5365.txt

% NCEP2_SIC=NCEP_SIC;
% NCEP2_SIC=NCEP2_SIC';
% % NCEP2_SIC(5476:5837)=[];
% NCEP2_SIC=reshape(NCEP2_SIC,4,365);
% NCEP2_SIC=mean(NCEP2_SIC);
% NCEP2_SIC=NCEP2_SIC';
% 
% E52_SIC=E5_SIC;
% E52_SIC=E52_SIC';
% %E52_SIC(5476:5837)=[];
% E52_SIC=reshape(E52_SIC,4,365);
% E52_SIC=mean(E52_SIC);
% E52_SIC=E52_SIC';

% load NCEP.txt
% 
figure 

% left = 0;
% right = left + 365;
% bottom = 0;
% top = bottom + 365;
% x = [left left right right];
% y = [bottom top top bottom];
% c = [0.20 0.11 0.13];
% fill(x, y, c,'FaceAlpha', 0.1)
% hold on

plot(SI125_NCEP365,'Color',[0.900, 0.50, 0.480],'LineWidth',2) %ncep
hold on

% plot(ASMR_SIC,'Color',[0.35, 0.40, 0.690],'LineWidth',2)
% hold on
% 
plot(SI125_ERA5365,'Color',[0, 0.75, 0.75],'LineWidth',2) %E5
hold on



title('Sea Ice Extent at 30°E,63°S','FontSize',25,'FontWeight','bold')
% legend('NCEP','ASMR','ERA5', 'FontSize',20,'FontWeight','bold',3)

xlim([304 365])
set(gca,'XTick',304:10:354)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% set(gca,'XTickLabel',{'jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',23,'FontWeight','bold')
set(gca,'XTickLabel',{'01-11-17','10-11-17','20-11-17','01-12-17','10-12-17','20-12-17'},'FontSize',23,'FontWeight','bold')

% set(gca,'YTick',-0.6:0.2:0.4)
% ylim([0 100])
ylabel('thickness[m]','FontSize',25,'FontWeight','bold')
% set(gca,'YTickLabel',{'0','10','20','30','40','50','60','70','80','90','100'},'FontSize',23,'FontWeight','bold')

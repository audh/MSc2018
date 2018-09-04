%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLAY ESIM RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANTARCTIC RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

% subplot(2,3,1)

plot(hs+hs_prec_bucket,'m','LineWidth',1.5)
hold on

plot(-hmi,'g','LineWidth',1.5) 
hold on
plot(-abs(hi_bio),'k','Linewidth',1)    %-abs(hi_bio) is plotting only the BAL values below the surface
hold on 

plot(-hmi-hi,'b','LineWidth',1.5)
hold on

title('Thickness and cover 1|ERA|8.0 W/m^2','FontSize',20,'FontWeight','bold')
legend('snow cover','meteoric ice','BAL','sea ice','FontSize',10,'FontWeight','bold')

xlim([1944 5837])
set(gca,'XTick',1944:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',20,'FontWeight','bold')
% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-1.0 0.3])
ylabel('thickness (m)','FontSize',20,'FontWeight','bold')
% set(gca,'YTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2'},'FontSize',10,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%% IRRADIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  subplot(2,3,2)
%  title('Irradiance')
%  
%  semilogy(I0,'m','LineWidth',1.5)
%  hold on
%  semilogy(IM,'g','LineWidth',1.5) 
%  hold on
%  semilogy(ISI,'b','LineWidth',1.5)
%  hold on
%  semilogy(ISI_bio,'k','LineWidth',1.5)
%  legend('surface I','interm layer I','sea ice I','BAL I','FontSize',10,'FontWeight','bold')
%  
%  xlim([0 5837])
% set(gca,'XTick',1:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% %  set(gca,'YTick',-0.6:0.2:0.4)
%  ylim([0.001 250])
%  ylabel('Irradiance (W m^-2)','FontSize',10,'FontWeight','bold')
% %  set(gca,'YTickLabel',{'0','50','100','150','200','250'},'FontSize',10,'FontWeight','bold')
%  set(gca,'Yscale','log')
% %%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,3)
% title('Temperature')
% 
% plot(Ta-273.15, 'y','LineWidth',1.2)
% hold on
% plot(T0-273.15,'r','LineWidth',1.2)
% hold on
% plot(Ts-273.15,'b','LineWidth',1.2) 
% hold on
% plot(Tsi-273.15,'g','LineWidth',1.2)
% hold on
% plot(Ti-273.15,'m','LineWidth',1.2)
% hold on
% plot(Tice_bio,'k','LineWidth',1.5)
% 
% legend('air T','surface T','snow T','interm layer T','sea ice T', 'BAL T', 'FontSize',5,'FontWeight','bold')
% 
% xlim([0 5837])
% set(gca,'XTick',1:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% 
% % set(gca,'YTick',-40:10:10)
% ylim([-30 20])
% ylabel('Temperature (deg celsius)','FontSize',10,'FontWeight','bold')
% % set(gca,'YTickLabel',{'-40','-20','0','20'},'FontSize',10,'FontWeight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BULK SALINITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% subplot(2,3,4)
% title('Bulk Salinity')
% 
% plot(Si,'b','LineWidth',1.2)
% hold on
% plot(Sice_bio,'k','LineWidth',1.2) 
% 
% legend('sea ice bulk S','BAL bulk S', 'FontSize',10,'FontWeight','bold')
% 
% xlim([0 5837])
% set(gca,'XTick',1:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% ylim([0 5])
% ylabel('Bulk Salinity (per mill)','FontSize',10,'FontWeight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BRINE SALINITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,5)
% title('Brine Salinity')
% 
% plot(Sbr_i,'b','LineWidth',1.2) 
% hold on
% plot(Sbr_bio,'k','LineWidth',1.2) 
% 
% legend('sea ice brine S','BAL brine S', 'FontSize',10,'FontWeight','bold')
% 
% xlim([0 5837])
% set(gca,'XTick',1:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% 
% ylabel('Brine Salinity (per mill)','FontSize',10,'FontWeight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BRINE VOLUME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,6)
% title('Brine Volume')
% 
% plot(Vbr_i,'b','LineWidth',1.2)
% hold on
% plot(Vbr_bio,'k','LineWidth',1.2)
% 
% 
% legend('sea ice brine V', 'BAL brine V', 'FontSize',10,'FontWeight','bold')
% 
% xlim([0 5837])
% set(gca,'XTick',1:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% 
% set(gca,'YTick',-0.6:0.2:0.4)
% ylim([0 0.2])
% ylabel('Relative Brine Volume','FontSize',10,'FontWeight','bold')
% set(gca,'YTickLabel',{'0.0','0.05','0.10','0.15','0.20'},'FontSize',10,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,6)
% figure
% title('BAL vs Time')
% 
% plot(hi_bio2,'b','LineWidth',1.2)
% hold on
% 
% legend('BAL thickness', 'FontSize',10,'FontWeight','bold')
% 
% xlim([0 5837])
%set(gca,'XTick',1:486.4167:5837)
%set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
% 
% set(gca,'YTick',-0.6:0.2:0.4)
% ylim([0 0.7])
% ylabel('thickness[m]','FontSize',10,'FontWeight','bold')
% set(gca,'YTickLabel',{'0.0','0.05','0.10','0.15','0.20'},'FontSize',10,'FontWeight','bold')

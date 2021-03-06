%% NCEP sensitivity test 

% N8a=hs+hs_prec_bucket;
% N8b=-hmi;
% N8c=-abs(hi_bio);
% N8d=-hmi-hi;

% N825a=hs+hs_prec_bucket;
% N825b=-hmi;
% N825c=-abs(hi_bio);
% N825d=-hmi-hi;

% N85a=hs+hs_prec_bucket;
% N85b=-hmi;
% N85c=-abs(hi_bio);
% N85d=-hmi-hi;

% N875a=hs+hs_prec_bucket;
% N875b=-hmi;
% N875c=-abs(hi_bio);
% N875d=-hmi-hi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(2,2,1)

plot(N8a,'m','LineWidth',1.5)
hold on

plot(N8b,'g','LineWidth',1.5) 
hold on
plot(N8c,'k','Linewidth',1)    %-abs(hi_bio) is plotting only the BAL values below the surface
hold on 

plot(N8d,'b','LineWidth',1.5)
hold on

title('Thickness and cover 1|NCEP|8.75 W/m^2','FontSize',15,'FontWeight','bold')
legend('snow cover','meteoric ice','BAL','sea ice','FontSize',3,'FontWeight','bold')

xlim([1944 5837])
set(gca,'XTick',1944:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',15,'FontWeight','bold')
% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-1.0 0.3])
ylabel('thickness (m)','FontSize',15,'FontWeight','bold')
% set(gca,'YTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2'},'FontSize',10,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%% IRRADIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure

subplot(2,2,2)

plot(N825a,'m','LineWidth',1.5)
hold on

plot(N825b,'g','LineWidth',1.5) 
hold on
plot(N825c,'k','Linewidth',1)    %-abs(hi_bio) is plotting only the BAL values below the surface
hold on 

plot(N825d,'b','LineWidth',1.5)
hold on

title('Thickness and cover 1|NCEP|8.75 W/m^2','FontSize',15,'FontWeight','bold')
legend('snow cover','meteoric ice','BAL','sea ice','FontSize',3,'FontWeight','bold')

xlim([1944 5837])
set(gca,'XTick',1944:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',15,'FontWeight','bold')
% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-1.0 0.3])
ylabel('thickness (m)','FontSize',15,'FontWeight','bold')
% set(gca,'YTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2'},'FontSize',10,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%% IRRADIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure

subplot(2,2,3)

plot(N85a,'m','LineWidth',1.5)
hold on

plot(N85b,'g','LineWidth',1.5) 
hold on
plot(N85c,'k','Linewidth',1)    %-abs(hi_bio) is plotting only the BAL values below the surface
hold on 

plot(N85d,'b','LineWidth',1.5)
hold on

title('Thickness and cover 1|NCEP|8.75 W/m^2','FontSize',15,'FontWeight','bold')
legend('snow cover','meteoric ice','BAL','sea ice','FontSize',3,'FontWeight','bold')

xlim([1944 5837])
set(gca,'XTick',1944:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',15,'FontWeight','bold')
% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-1.0 0.3])
ylabel('thickness (m)','FontSize',15,'FontWeight','bold')
% set(gca,'YTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2'},'FontSize',10,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%% IRRADIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure

subplot(2,2,4)

plot(N875a,'m','LineWidth',1.5)
hold on

plot(N875b,'g','LineWidth',1.5) 
hold on
plot(N875c,'k','Linewidth',1)    %-abs(hi_bio) is plotting only the BAL values below the surface
hold on 

plot(N875d,'b','LineWidth',1.5)
hold on

title('Thickness and cover 1|NCEP|8.75 W/m^2','FontSize',15,'FontWeight','bold')
legend('snow cover','meteoric ice','BAL','sea ice','FontSize',3)

xlim([1944 5837])
set(gca,'XTick',1944:486.4167:5837)
% set(gca,'XTickLabel',{'jan-17','feb-17','mar-17','apr-17','may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',10,'FontWeight','bold')
set(gca,'XTickLabel',{'may-17','jun-17','jul-17','aug-17','sep-17','oct-17','nov-17','dec-17'},'FontSize',15,'FontWeight','bold')
% set(gca,'YTick',-0.6:0.2:0.4)
ylim([-1.0 0.3])
ylabel('thickness (m)','FontSize',15,'FontWeight','bold')
% set(gca,'YTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2'},'FontSize',10,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%% IRRADIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

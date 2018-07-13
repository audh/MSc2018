%% MATLAB FINAL
% reading files
file = dir('/home/student/ADHRIE001/data/A2017*.nc');
log_index = [file.isdir]; %creates a logical index
filename = {file(~log_index).name}; %gets file name
%%
% for looooooop

for i = 1:4 %index of loop
    ncfile = char(filename(i)); %accesses the ith file name and converts to char array
    lon = ncread(ncfile,'lon');
    lat = ncread(ncfile,'lat');
    chlo = ncread(ncfile,'chlor_a');
    lon_index = find(lon>=10 & lon<=40);
    lat_index = find(lat>=-50 & lat<=-30);      %extracting indecies of non-zero elements
    chlo_data(i,:,:) = chlo(lon_index,lat_index);      
end
% 2. nanmean
mean_chlo = squeeze(nanmean(chlo_data,1)); %nanmean disregards the NaN values, squeeze removes the singleton dimension

% 3. climatology data

clim_data = 'A20023352016366.L3m_MC_CHL_chlor_a_9km.nc';
%ncdisp(clim_data)
lon_c = ncread(clim_data,'lon');
lat_c = ncread(clim_data,'lat');
chlo_c = ncread(clim_data,'chlor_a');
lon_index_c = find(lon_c>=10 & lon_c<=40);
lat_index_c = find(lat_c>=-50 & lat_c<=-30);      %extracting indecies of non-zero elements
chlo_data_c = chlo_c(lon_index_c,lat_index_c);

% 4. anomaly
anomaly = mean_chlo - chlo_data_c ;

% 5. Plot
% plot 1: Jan comp.

figure
subplot(1,2,1)       % add first plot in 1 x 2 grid
m_proj('robinson',...        %info about type of projection
    'lon',[10 40],'lat',[-50  -30])
%m_grid
m_coast('color',[0 0 0])        %colour is just the coast line only, patch (used in sa map) is a solid continent
m_grid                          % grid is added to figure, no specific details (as in sa map)
hold on                        %hold ON holds the current plot and all axis properties so that subsequent graphing commands add to the existing graph.
m_pcolor(lon(lon_index),lat(lat_index),mean_chlo')          %assigns colour scheme        
shading flat                    % shading FLAT sets the shading of the current graph to flat.
colorbar                        %displays colour bar
caxis([0 5])                  %scales axis
title('Composite','Fontweight','bold')
%print -dpng -r300 ADHRIE001_JanC1.png;

% figure 2: anomaly plot

subplot(1,2,2)       % add second plot in 2 x 1 grid
m_proj('robinson',...        %info about type of projection
    'lon',[10 40],'lat',[-50  -30])
%m_grid
m_coast('color',[0 0 0])        %colour is just the coast line only, patch (used in sa map) is a solid continent
m_grid                          % grid is added to figure, no specific details (as in sa map)
hold on                         %hold ON holds the current plot and all axis properties so that subsequent graphing commands add to the existing graph.
m_pcolor(lon(lon_index_c),lat(lat_index_c),anomaly')          %assigns colour scheme        
shading flat                    % shading FLAT sets the shading of the current graph to flat.
colorbar                        %displays colour bar
caxis([-1 1])                  %scales c axis
colormap(bluewhitered)
title('Anomaly','Fontweight','bold')
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'Chlorophyll Concentration for January 2017','HorizontalAlignment','center','VerticalAlignment', 'top','Fontweight','bold')
print -dpng -r300 ADHRIE001_JanC2.png;

%% end %%

% extract data from a point in .nc file %

% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/air.2m.gauss.2017.nc' ;
% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/dswrf.sfc.gauss.2017.nc' ;
% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/prate.sfc.gauss.2017.nc' ;
% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/shum.2m.gauss.2017.nc' ;
% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/tcdc.eatm.gauss.2017.nc' ;
% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/uwnd.10m.gauss.2017.nc' ;
% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/vwnd.10m.gauss.2017.nc' ;

% file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/skt.sfc.gauss.2017.nc' ;
file='/home/riesna/MSc2018/ESIM_AntTest/ncfiles/rhum.sig995.2017.nc' ;
ncdisp(file)

longitude=ncread(file,'lon');
latitude=ncread(file,'lat');
time=ncread(file,'time');
rhum=ncread(file,'rhum');
% level=ncread(file,'level');
a=[62,61,60];
% a=[81,80,79];
 basename='rh';

for i = 1:3 %index of loop
%    skt= ncread(file,'skt',[17 a(i) 1 1],[1 1 1 inf],[1 1 1 1]);
   rhum= ncread(file,'rhum',[13 a(i) 1],[1 1 inf],[1 1 1]);
   rh= reshape(squeeze(rhum(:,1,:,:)),1460,1);
   filename=[basename,num2str(i),'.txt'];
   dlmwrite(filename,rh);   
end


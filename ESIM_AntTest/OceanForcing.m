%% Ocean forcing file sorting %%

TEMP=textread('/home/riesna/MSc2018/ESIM_AntTest/temp1.txt'); %#ok<DTXTRD>
SAL=textread('/home/riesna/MSc2018/ESIM_AntTest/sal1.txt'); %#ok<DTXTRD>
day=736696:1:737060; % 01-01-17 : 31:12:17 (365d)		This is specific for Kobbefjord test case study site. To be adjusted by the user as needed.
day=day';
day=datestr(day,31);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CMEMS_forcing_num=[TEMP SAL];

CMEMS_forcing_str=num2str(CMEMS_forcing_num);

time=num2str(day);

BFM_input_ocean_ANT=[time CMEMS_forcing_str];

dlmwrite('BFM_input_ocean_ANT.dat', BFM_input_ocean_ANT, 'delimiter', '%s',  'precision', 6)

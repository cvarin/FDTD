% Copyright (c) 2012, Ilker R Capoglu 
% This script reads an Angora line file and plots the data.

fclose('all');
clear all;

linefilename = '/my_path/LineFile_Ey_Y_0_1.aln';

linefile = fopen(linefilename);

version_major = fread(linefile ,1,'int');
version_minor = fread(linefile ,1,'int');
version_rev = fread(linefile ,1,'int');
dt = fread(linefile ,1,'double');
initial_time_value = fread(linefile ,1,'double');
total_length = fread(linefile ,1,'int');
length_time = fread(linefile ,1,'int');
PML = fread(linefile ,1,'int');

plot_length = total_length-2*PML;

fullarray = zeros(1,total_length);
clippedarray = zeros(1,plot_length);
linefig = figure;
lineaxes = axes;
lineplot=plot(lineaxes,(1:plot_length),clippedarray,'YDataSource','clippedarray');
maxvalue = 1;
minvalue = -1;
axis([1 plot_length minvalue maxvalue]);

pause_time = 0.01;

for n=1:length_time
    fullarray = fread(linefile,total_length,'double');
    clippedarray = fullarray(PML+1:total_length-PML);
    refreshdata(lineplot);
    drawnow;
    pause(pause_time);
    title(['Time step : ',num2str(n)]);
end

fclose(linefile);
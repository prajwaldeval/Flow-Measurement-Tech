clc;
clear all;

%%Load PIV files
load('0mean.mat');
PIVm = Vloc;
load('0rms.mat');
PIVrms = Vloc;
clear Vloc;
load('yPIV.mat')

%%Load HWA file
load('hwavel_00.mat');
HWA = mean_velocity_meas;
yFMT = -40:4:40;
yFMT(end) = 39;

%%plotting
plot(PIVm,yloc,'-', 'color','cyan');
hold on
plot(PIVrms,yloc, '-', 'color','red');
plot(HWA(2:end),yFMT(1:end-2),'color','black');
clc;
clear all;

%%Define Location of files
path = 'Data/Drive/';
subfolder = '15_optimal';
meanrms = strcat('/',subfolder,'_stats_10_images');

%%Same for all
instant = '/B00001.dat';
mean = '/B00001.dat';
rms = '/B00002.dat';

%%Reading
fullpath = [path,subfolder,meanrms,mean];
fid = fopen(fullpath);
data = textscan(fid, '%f %f %f %f', 'headerlines', 3);
x = data{1};
y = data{2};
u = data{3};
v = data{4};
V = sqrt(u.*u + v.*v);
clear data;

    
% scatter(x,y,[],V,'filled')
% colorbar
% colormap('gray')


%%Extract c/t=1.2
loc = x==120.000601000000;
xloc = x(loc);
yloc = y(loc);
uloc = u(loc);
vloc = v(loc);
Vloc = V(loc);

plot(Vloc,yloc);

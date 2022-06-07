clc;
clear all;

% Reading TIF images
im_cal = imread('Data/B00001_cal.tif');

im = imread('Data/B00010.tif');
[size_y,size_x] = size(im);
height = size_y/2;
im1 = im(1:height,:);
im2 = im(height + 1:end,:);
clear im;

% Background computation
im_min = ((im1+im2) - abs(im1-im2))/2;
im_ave = (im1 + im2)/2;
im_sub = im1 - im2;

[xmax, ymax] = size(im1);

% Window size 
wsize = [32, 32];
w_width = wsize(1);
w_height = wsize(2);

% Center points grid 
xmin = w_width/2;
ymin = w_height/2;
xgrid = 200:w_width/2:864;
ygrid = 200:w_height/2:1696;

% Number of windows in total
w_xcount = length(xgrid);
w_ycount = length(ygrid);

% Range for search windows
x_disp_max = w_width/2;
y_disp_max = w_height/2;

xpeakl = 0;
ypeakl = 0;

test_im1(w_width, w_height) = 0; 
test_im2(w_width+2*x_disp_max, w_height + 2 * y_disp_max) = 0; 
dpx(w_xcount, w_ycount) = 0; 
dpy(w_xcount, w_ycount) = 0; 

xpeak1 = 0; 
ypeak1 = 0;

% i, j are for the windows
% test_i and test_j are for the test window to be
% extracted from image A
for i=1:(w_xcount)
    for j=1:(w_ycount)
        max_correlation = 0;
        test_xmin = xgrid(i) - w_width/ 2;
        test_xmax = xgrid(i) + w_width/ 2;
        test_ymin = ygrid(j) - w_height/2;
        test_ymax = ygrid(j) + w_height/2;
        x_disp = 0;
        y_disp = 0;
        test_im1 = im1(test_xmin:test_xmax, test_ymin:test_ymax);
        test_im2 = im2(test_xmin:test_xmax, test_ymin:test_ymax);
        correlation = normxcorr2(test_im1,test_im2);
        [xpeak, ypeak] = find (correlation == max(correlation(:)));
        
        % Re-scaling
        xpeakl = test_xmin + xpeak - wsize(1)/2 - x_disp_max;
        ypeakl = test_ymin + ypeak - wsize(2)/2 - y_disp_max;
        dpx(i,j) = xpeakl - xgrid(i);
        dpy(i,j) = ypeakl - ygrid(j);
        
    end
end

%Vector Display 
quiver(dpx, - dpy)



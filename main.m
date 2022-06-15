clc;
clear all;

%some fax
pxtocm = 94.9; %px/cm from cal image
dt = 70e-6; %time between the two images


%set vars
windowsize = 24;
search_size = 2;%as proportion of window size
stdev_threshold= 1;

% Reading TIF images
im_cal = imread('Data/B00001_cal.tif');

im = imread('Data/B00010.tif');
[size_y,size_x] = size(im);
size_y = size_y/2; %pictures are stacked together so height is actually half
im1 = im(1:size_y,:);
im2 = im(size_y + 1:end,:);
clear im;



% Window size 
wsize = [windowsize, windowsize];
w_height = wsize(1);
w_width = wsize(2);

% Search Window Size
x_search = round(w_width*search_size);
y_search = round(w_height*search_size);

% Center points grid
xmin = w_width/2;
ymin = w_height/2;
xgrid = w_width:w_width:size_x;  %center point of each window in x
ygrid = w_height:w_height:size_y;%center point of each window in y

% Number of windows in total
w_xcount = length(xgrid);
w_ycount = length(ygrid);



%%Declaring vars for loop
window(w_height, w_width) = 0; %initialise subwindow var
reference( w_height + 2 * y_search, w_width+2*x_search) = 0; %initialise search area var
dpx(w_ycount, w_xcount) = 0; %initialise vars for velocity comp in x
dpy(w_ycount, w_xcount) = 0; %initialise vars for velocity comp in y


im1bg = im1;
im2bg = im2;
for i=1:(w_xcount)
    for j=1:(w_ycount)
        max_correlation = 0;
        test_xmin = xgrid(i) - w_width/ 2 + 1;
%         if test_xmin < 1
%             test_xmit = 1;
%         end
        
        test_xmax = xgrid(i) + w_width/ 2;
        if test_xmax > size_x
            test_xmax = size_x;
        end
       
        test_ymin = ygrid(j) - w_height/2 + 1;
%         if test_ymin < 1
%             test_ymin = 1;
%         end
        
        test_ymax = ygrid(j) + w_height/2;
        if test_ymax > size_y
            test_ymax = size_y;
        end
        
        x_disp = 0;
        y_disp = 0;
        
        
        
        window1 = im1(test_ymin:test_ymax, test_xmin:test_xmax); %32x32 test window
        window2 = im2(test_ymin:test_ymax, test_xmin:test_xmax); %32x32 test window (only for bg subtraction purposes)
        
        bg1 = min(window1,[],'all');
        bg2 = min(window2,[],'all');
        bg = min([bg1,bg2]);
        bgr1 = window1 - bg; %window 1 with bg removed
        bgr2 = window2 - bg; %window 2 with bg removed

        im1(test_ymin:test_ymax, test_xmin:test_xmax) = bgr1;
        im2(test_ymin:test_ymax, test_xmin:test_xmax) = bgr2;

    end
end

%looping through each 32x32 window to find its most likely velocity
%component

for i=1:(w_xcount)
    for j=1:(w_ycount)
        max_correlation = 0;
        test_xmin = xgrid(i) - w_width/ 2 + 1;
%         if test_xmin < 1
%             test_xmit = 1;
%         end
        
        test_xmax = xgrid(i) + w_width/ 2;
        if test_xmax > size_x
            test_xmax = size_x;
        end
       
        test_ymin = ygrid(j) - w_height/2 + 1;
%         if test_ymin < 1
%             test_ymin = 1;
%         end
        
        test_ymax = ygrid(j) + w_height/2;
        if test_ymax > size_y
            test_ymax = size_y;
        end
        
        x_disp = 0;
        y_disp = 0;
        
        
        
        test = im1(test_ymin:test_ymax, test_xmin:test_xmax); %32x32 test window
        
        
        ref_xmin = test_xmin - x_search + 1;
        if ref_xmin < 1
            ref_xmin = 1;
        end
        
        ref_xmax = test_xmax + x_search;
        if ref_xmax > size_x
            ref_xmax = size_x;
        end
        
        ref_ymin = test_ymin - y_search + 1;
        if ref_ymin < 1
            ref_ymin = 1;
        end
        
        ref_ymax = test_ymax + y_search;
        if ref_ymax > size_y
            ref_ymax = size_y;
        end
        
        reference = im2(ref_ymin:ref_ymax,ref_xmin:ref_xmax);
        
        if sum(test,'all') == 0
%             disp('Assume 0');
            dpx(j,i) = 0;
            dpy(j,i) = 0;
        else
%           disp('Calculating xcorr');
          correlation = normxcorr2(test,reference);
          [maxcorr_y, maxcorr_x] = find(abs(correlation) == max(abs(correlation(:))));
          dpx(j,i) = ref_xmin - w_width/2 + maxcorr_x - xgrid(i);
          dpy(j,i) = ref_ymin - w_height/2 + maxcorr_y - ygrid(j);
        end
        
        
        % Re-scaling
%         xpeakl = test_xmin + xpeak - wsize(1)/2 - x_search;
%         ypeakl = test_ymin + ypeak - wsize(2)/2 - y_search;
%         dpx(i,j) = xpeakl - xgrid(i);
%         dpy(i,j) = ypeakl - ygrid(j);
        
    end
end
stdev_threshold = 1;
vx = fliplr(0.01*((dpx./dt) /pxtocm));
vy = fliplr(0.01*((dpy./dt) /pxtocm));
v = sqrt(vx.*vx + vy.*vy);
stdev = std(v,0,'all');
meanv = mean(mean(v));

vx((v-meanv)/stdev > stdev_threshold) = NaN;
vy((v-meanv)/stdev > stdev_threshold) = NaN;

quiver(vx,vy);

%Vector Display 
% quiver(dpx, - dpy)



clc;
clear all;


%set vars
windowsize = 32;
search_size = 0.75;%as proportion of window size

% Reading TIF images
im_cal = imread('Data/B00001_cal.tif');

im = imread('Data/B00010.tif');
[size_y,size_x] = size(im);
size_y= size_y/2; %pictures are stacked together so height is actually half
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
xgrid = w_width/2:w_width:size_x;  %center point of each window in x
ygrid = w_height/2:w_height:size_y;%center point of each window in y

% Number of windows in total
w_xcount = length(xgrid);
w_ycount = length(ygrid);

xpeakl = 0;
ypeakl = 0;

%%Declaring vars for loop
window(w_width, w_height) = 0; %initialise subwindow var
reference(w_width+2*x_search, w_height + 2 * y_search) = 0; %initialise search area var
dpx(w_xcount, w_ycount) = 0; %initialise vars for velocity comp in x
dpy(w_xcount, w_ycount) = 0; %initialise vars for velocity comp in y


im1bg = im1;
im2bg = im2;


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
        
        
        
        window1 = im1(test_ymin:test_ymax, test_xmin:test_xmax); %32x32 test window
        window2 = im2(test_ymin:test_ymax, test_xmin:test_xmax); %32x32 test window (only for bg subtraction purposes)
        
        bg1 = min(window1,[],'all');
        bg2 = min(window2,[],'all');
        bg = min([bg1,bg2]);
        bgr1 = window1 - bg; %window 1 with bg removed
        bgr2 = window2 - bg; %window 2 with bg removed

        im1(test_ymin:test_ymax, test_xmin:test_xmax) = bgr1;
        im2(test_ymin:test_ymax, test_xmin:test_xmax) = bgr2;
        
        
        ref_xmin = test_xmin - x_search/2 + 1;
        if ref_xmin < 1
            ref_xmin = 1;
        end
        
        ref_xmax = test_xmax + x_search/2;
        if ref_xmax > size_x
            ref_xmax = size_x;
        end
        
        ref_ymin = test_ymin - y_search/2 + 1;
        if ref_ymin < 1
            ref_ymin = 1;
        end
        
        ref_ymax = test_ymax + y_search/2;
        if ref_ymax > size_y
            ref_ymax = size_y;
        end
        
        reference = im2(ref_ymin:ref_ymax,ref_xmin:ref_xmax);
        
        if sum(bgr1,'all') == 0
            dpx(i,j) = 0;
            dpy(i,j) = 0;
        else
          correlation = normxcorr2(bgr1,reference);
        end
        
  
%         correlation = normxcorr2(bgr1,reference);
%         [xpeak, ypeak] = find (correlation == max(correlation(:)));
        
        % Re-scaling
%         xpeakl = test_xmin + xpeak - wsize(1)/2 - x_search;
%         ypeakl = test_ymin + ypeak - wsize(2)/2 - y_search;
%         dpx(i,j) = xpeakl - xgrid(i);
%         dpy(i,j) = ypeakl - ygrid(j);
        
    end
end

%Vector Display 
% quiver(dpx, - dpy)



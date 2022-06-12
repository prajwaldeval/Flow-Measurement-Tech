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

% Background computation, different methodologies - none of these really
% make sense
im_min = ((im1+im2) - abs(im1-im2))/2;
im_ave = (im1 + im2)/2;
im_sub = im1 - im2;


% Window size 
wsize = [windowsize, windowsize];
w_height = wsize(1);
w_width = wsize(2);

% Search Window Size
x_search = w_width*search_size;
y_search = w_height*search_size;

% Center points grid 
xmin = w_width/2;
ymin = w_height/2;
xgrid = x_search:w_width/2:size_x-x_search; %starts & ends x&y grid such that
ygrid = y_search:w_height/2:size_y-y_search;%there is always a search
                                            %window sized margin outside


                                            
% % % % % CHANGES FIXED UPTO HERE
                                            
% Number of windows in total
w_xcount = length(xgrid);
w_ycount = length(ygrid);



xpeakl = 0;
ypeakl = 0;

test_im1(w_width, w_height) = 0; 
test_im2(w_width+2*x_search, w_height + 2 * y_search) = 0; 
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
        xpeakl = test_xmin + xpeak - wsize(1)/2 - x_search;
        ypeakl = test_ymin + ypeak - wsize(2)/2 - y_search;
        dpx(i,j) = xpeakl - xgrid(i);
        dpy(i,j) = ypeakl - ygrid(j);
        
    end
end

%Vector Display 
quiver(dpx, - dpy)



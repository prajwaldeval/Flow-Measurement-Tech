clc;
clear all;

im_cal = imread('Data/B00001_cal.tif');

im = imread('Data/B00010.tif');
[size_y,size_x] = size(im);
height = size_y/2;
im1 = im(1:height,:);
im2 = im(height + 1:end,:);
clear im;
im_min = ((im1+im2) - abs(im1-im2))/2;
im_ave = (im1 + im2)/2;
im_sub = im1 - im2;
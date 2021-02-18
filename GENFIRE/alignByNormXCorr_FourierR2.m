%%  alignByNormXCorr %%

%% Align projections using normalized cross correlation
%%inputs:
%%  ref_img           - reference image; corresponds to the input projection
%%  img               - image to compare; corresponds to the calculated backprojection
%%  smooth_factor     - smoothing parameter for smooth3D. The backprojection is smoothed and thresholded to form a template
%%  threshhold_factor - value between 0 and 1 representing the fraction of maximum intensity in backprojection to use as a threshhold

%%outputs:
%%  XC - the maximum value of cross correlation
%%  new_x_center - the x (dimension 1) location of the best template match in ref_img
%%  new_y_center - the y (dimension 2) location of the best template match in ref_img

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function [XC, new_x_center, new_y_center] = alignByNormXCorr_FourierR2(ref_img,img,smooth_factor,threshhold_factor)

if nargin<3
    smooth_factor = 0.25;
    threshhold_factor = 0.10;
end

[dimx, dimy] = size(ref_img);
ncx = round((dimx+1)/2);
ncy = round((dimy+1)/2);

if nargin<3
    window_half_size = floor(dimx/10);
end

bestXC = 0;
new_x_center = ncx;
new_y_center = ncy;

simg = smooth3D(img,smooth_factor);
threshhold = max(simg(:)) .* threshhold_factor;
m = smooth3D(img,smooth_factor)>threshhold;
[sx,sy] = find(m~=0);
template = img(min(sx(:)):max(sx(:)),min(sy(:)):max(sy(:)));
template=img;
xc = normxcorr2(template,ref_img);
%XC = max(xc(:));

[xpeak, ypeak]   = find(xc==max(xc(:)));
x_offSet = xpeak-size(template,1);
y_offSet = ypeak-size(template,2);
new_x_center = x_offSet + round((size(template,1)+1)/2);
new_y_center = y_offSet + round((size(template,2)+1)/2);

shiftX = (ncx - new_x_center);
shiftY = (ncy - new_y_center);
    
ref_imgShifted = circshift(ref_img,[shiftX, shiftY]);

FF_ref = My_FFTN(ref_imgShifted);
FF_img = My_FFTN(img);
XC = sum( (FF_ref(:)-FF_img(:)).^2);

end
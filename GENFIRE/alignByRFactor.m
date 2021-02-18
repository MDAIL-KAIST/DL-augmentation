%%  alignByRFactor %%

%% Align projections using real space R factor
%%inputs:
%%  ref_img           - reference image; corresponds to the input projection
%%  img               - image to compare; corresponds to the calculated backprojection
%%  window_half_size  - R factor will be calculated for translations of +- window_half_size

%%outputs:
%%  R            - the minimum R factor found
%%  new_x_center - the x (dimension 1) location of the best match in ref_img
%%  new_y_center - the y (dimension 2) location of the best match in ref_img

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function [R, new_x_center, new_y_center] = alignByRFactor(ref_img,img,window_half_size)

if all(size(ref_img)~=size(img))
    error('Images must be same size')
end

[dimx, dimy] = size(ref_img);
ncx = round((dimx+1)/2);
ncy = round((dimy+1)/2);

if nargin<3
    window_half_size = floor(dimx/10);
end

bestR = 1e30;
new_x_center = ncx;
new_y_center = ncy;
norm_factor = sum(abs(ref_img(:)));
for sx = -window_half_size:window_half_size
    for sy = -window_half_size:window_half_size
        R = sum(sum(abs(ref_img-circshift(img,[sx, sy])))) ./ norm_factor;
        if R < bestR
           bestR = R;
           new_x_center = ncx - sx;
           new_y_center = ncy - sy;
        end
    end
end

end
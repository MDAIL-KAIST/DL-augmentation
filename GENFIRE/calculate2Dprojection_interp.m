%%  calculate2Dprojection_interp %%

%%calculates 2D projection from 3D object using linear interpolation of
%%central slice in Fourier space
%%inputs:
%%  modelK - 3D Fourier space of object to compute projection of
%%  phi - rotation angle

%%outputs:
%%  projection - result

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.



function [projection, kSlice] = calculate2Dprojection_interp(modelK,phi)

%get dimensions and centers
[dimx, dimy] = size(modelK);

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 

X = (1:dimx) - ncx;
Y = zeros(size(X));
%calculate rotation matrix
c = cosd(phi);
s = sind(phi);
R = [c, -s;
     s, c];

  
% [ky, kx] = meshgrid((1:dimy) - ncy, (1:dimx) - ncx);

%rotate coordinates
% rotkCoords = R*[X;Y];
rotkCoords = R*[X;Y];
rotKx = rotkCoords(1,:);
rotKy = rotkCoords(2,:);

%reshape for interpolation
rotKx = reshape(rotKx,size(X));
rotKy = reshape(rotKy,size(Y));

%calculate points on central slice
% kSlice = interp2(ky,kx,modelK,rotKy,rotKx,'linear');
% kSlice = interp2(ky,kx,modelK,rotKy,rotKx,'cubic');
kSlice = splinterp2(modelK,rotKy+ncx,rotKx+ncx);

%remove any nan from interpolation
kSlice(isnan(kSlice))=0;

%take IFFT to obtain projection
projection = kSlice;
% projection = real(my_ifft(kSlice(:,:)));
end
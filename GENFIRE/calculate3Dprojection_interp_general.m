%%  calculate3Dprojection_interp %%

%%calculates 2D projection from 3D object using linear interpolation of
%%central slice in Fourier space
%%inputs:
%%  modelK - 3D Fourier space of object to compute projection of
%%  phi,theta,psi - Euler angles of desired projection

%%outputs:
%%  projection - result

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.



function [projection, kSlice] = calculate3Dprojection_interp_general(modelK,phi,theta,psi,vec1,vec2,vec3)

%get dimensions and centers
[dimx, dimy, dimz] = size(modelK);

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
ncz = round((dimz+1)/2);

[Y, X, Z] = meshgrid( ((1:dimy) - ncy) / dimy, ((1:dimx) - ncx) / dimx, 0);

%calculate rotation matrix
% R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
%       -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
%       sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];
R1 = MatrixQuaternionRot(vec1,phi);
R2 = MatrixQuaternionRot(vec2,theta);
R3 = MatrixQuaternionRot(vec3,psi);
% R1 = [cosd(phi) -sind(phi) 0;
%       sind(phi) cosd(phi) 0;
%       0          0         1];    
%     
% R2 = [cosd(theta) 0  sind(theta);
%         0         1      0;
%       -sind(theta) 0 cosd(theta)];
% 
% R3 = [1   0   0;
%        0  cosd(psi) -sind(psi);
%       0   sind(psi) cosd(psi)];    
    
 R =(R1*R2*R3)';

[ky kx kz ] = meshgrid(((1:dimy) - ncy) / dimy, ((1:dimx) - ncx) /  dimx, ((1:dimz) - ncz) /dimz);

%rotate coordinates
rotkCoords = R'*[X(:)';Y(:)';Z(:)'];
rotKx = rotkCoords(1,:);
rotKy = rotkCoords(2,:);
rotKz = rotkCoords(3,:);

%reshape for interpolation
rotKx = reshape(rotKx,size(X));
rotKy = reshape(rotKy,size(Y));
rotKz = reshape(rotKz,size(Z));

%calculate points on central slice
% kSlice = interp3(ky,kx,kz,modelK,rotKy,rotKx,rotKz,'linear');
% kSlice = interp3(ky,kx,kz,modelK,rotKy,rotKx,rotKz,'cubic');
%kSlice = splinterp3(modelK,(rotKy+ncy/dimy)*dimy,(rotKx+ncx/dimx)*dimx,(rotKz+ncz/dimz)*dimz);
kSlice = interp3(ky,kx,kz,modelK,rotKy,rotKx,rotKz);
%remove any nan from interpolation
kSlice(isnan(kSlice))=0;

%take IFFT to obtain projection
projection = real(my_ifft(kSlice(:,:)));
end
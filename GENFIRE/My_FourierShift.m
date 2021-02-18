% My_FourierShift

% Author: Yongsoo Yang, UCLA Physics & Astronomy
%         yongsoo.ysyang@gmail.com

% This code runs sub-pixel circular shift based on Fourier factoring

function img2 = My_FourierShift(img,dy,dx)

ny = size(img,1); nx = size(img,2);
[X, Y] = meshgrid(-ceil((nx-1)/2):floor((nx-1)/2),-ceil((ny-1)/2):floor((ny-1)/2));
F = fftshift(ifftn(ifftshift(img)));
Pfactor = exp(2*pi*1i*(dx*X/nx + dy*Y/ny));
img2 = fftshift(fftn(ifftshift(F.*Pfactor)));

end
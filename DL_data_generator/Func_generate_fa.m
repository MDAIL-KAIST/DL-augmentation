% J.Lee, KAIST, 2020
% input : AtomNumbers, volsize
% output : fa (electron scattering factor in 3D fourier space, volsize^3)

function fa = Func_generate_fa(AtomNumbers,volsize,Res)

if length(volsize) == 3
    sizeX = [volsize(1) volsize(2) volsize(3)];
elseif length(volsize) == 1
    sizeX = [volsize volsize volsize];
else
    error('volsize should be either length 3 or length 1!')
end

fa = zeros(sizeX(1),sizeX(2),sizeX(3),length(AtomNumbers));

finalvol_summed = zeros(sizeX);

kx = 1:size(finalvol_summed,1);
ky = 1:size(finalvol_summed,2);
kz = 1:size(finalvol_summed,3);

MultF_X = 1/(length(kx)*Res);
MultF_Y = 1/(length(ky)*Res);
MultF_Z = 1/(length(kz)*Res);

CentPos = round((size(finalvol_summed)+1)/2);
[KX, KY, KZ] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y,(kz-CentPos(3))*MultF_Z);
q2 = KX.^2 + KY.^2 + KZ.^2;
clear KX KY KZ

for j=1:length(AtomNumbers)
    fa(:,:,:,j) = reshape(fatom_vector( sqrt(q2),AtomNumbers(j)), sizeX);
end

end
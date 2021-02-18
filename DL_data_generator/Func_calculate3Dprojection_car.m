function Proj_car = Func_calculate3Dprojection_car(Vol,phi,theta,psi,vec1,vec2,vec3)
% J.Lee, KAIST, 2020

%get dimensions and centers
[dimx, dimy, dimz] = size(Vol);

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
ncz = round((dimz+1)/2);

%calculate rotation matrix
R1 = MatrixQuaternionRot(vec1,phi);
R2 = MatrixQuaternionRot(vec2,theta);
R3 = MatrixQuaternionRot(vec3,psi);
  
R =(R1*R2*R3)';

[Ry Rx Rz ] = meshgrid(((1:dimy) - ncy) / dimy, ((1:dimx) - ncx) /  dimx, ((1:dimz) - ncz) /dimz);

%rotate coordinates
rotRCoords = R'*[Rx(:)';Ry(:)';Rz(:)'];
rotRx = rotRCoords(1,:);
rotRy = rotRCoords(2,:);
rotRz = rotRCoords(3,:);

%reshape for interpolation
rotRx = reshape(rotRx,size(Rx));
rotRy = reshape(rotRy,size(Ry));
rotRz = reshape(rotRz,size(Rz));

%rotation & summation z component
RSlice = interp3(Ry,Rx,Rz,Vol,rotRy,rotRx,rotRz);
RSlice(isnan(RSlice))=0;
RSlice(RSlice<0)=0;
Proj_car = sum(RSlice,3);
%Proj_car = permute(sum(RSlice,1),[2 3 1]);
end
function Proj = calculate3Dprojection_RealSpaceinterp(Vol, phi, theta, psi)
    
%     vector1 = [0 0 1];
%     rotmat1 = MatrixQuaternionRot(vector1,phi);    
% 
%     vector2 = [0 1 0];
%     rotmat2 = MatrixQuaternionRot(vector2,theta);
% 
%     vector3 = [0 0 1];
%     rotmat3 = MatrixQuaternionRot(vector3,psi);
% 
%     R =  (rotmat1*rotmat2*rotmat3);
R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,  -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi) ,  sind(theta)*cosd(phi)  ;
      cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,  sind(theta)*sind(phi)   ;
               -cosd(psi)*sind(theta)    ,                          sind(psi)*sind(theta)       ,              cosd(theta)]; 

    VolSize = size(Vol);
    Xarr = -ceil((VolSize(1)-1)/2):floor((VolSize(1)-1)/2);
    Yarr = -ceil((VolSize(2)-1)/2):floor((VolSize(2)-1)/2);
    Zarr = -ceil((VolSize(3)-1)/2):floor((VolSize(3)-1)/2);

    [Y, X, Z] = meshgrid(Yarr, Xarr, Zarr);
    
    Shape = size(X);

    posvec = [X(:), Y(:), Z(:)]';
    rot_posvec = R*posvec;

    rotX = reshape(rot_posvec(1,:),Shape);
    rotY = reshape(rot_posvec(2,:),Shape);
    rotZ = reshape(rot_posvec(3,:),Shape);

    clear posvec;
    clear rot_posvec;

    cut_size = 64;

    cutnum = floor(Shape(3)/cut_size);
    cutrem = mod(Shape(3),cut_size);

    ROTvol = zeros(Shape);
    
    for i=1:cutnum
        cutindar = ((i-1)*cut_size+1):i*cut_size;
        ROTtemp = interp3(Y,X,Z,Vol,rotY(:,:,cutindar),rotX(:,:,cutindar),rotZ(:,:,cutindar),'cubic',0);
        ROTvol(:,:,cutindar) = ROTtemp;    
    end

    if cutrem~=0
        cutindar = (cutnum*cut_size+1):size(rotX,3);
        ROTtemp = interp3(Y,X,Z,Vol,rotY(:,:,cutindar),rotX(:,:,cutindar),rotZ(:,:,cutindar),'cubic',0);
        ROTvol(:,:,cutindar) = ROTtemp;
    end
        
    Proj = squeeze(sum(ROTvol,3));

end



% 
% function dd = MatrixQuaternionRot(vector,theta)
% 
%   theta = theta*pi/180;
%   vector = vector/sqrt(dot(vector,vector));
%   w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
%   RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
%         2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
%         2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];
% 
%   dd = RotM;
% end
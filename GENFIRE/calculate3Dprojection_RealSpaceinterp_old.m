function Proj = calculate3Dprojection_RealSpaceinterp(Vol, phi, theta, psi)

    vector1 = [0 0 1];
    rotmat1 = MatrixQuaternionRot(vector1,phi);    

    vector2 = [0 1 0];
    rotmat2 = MatrixQuaternionRot(vector2,theta);

    vector3 = [0 0 1];
    rotmat3 = MatrixQuaternionRot(vector3,psi);

    rotMATs =  (rotmat1*rotmat2*rotmat3);
    
    Xarr = My_find_symm_indarr(size(Vol,1));
    Yarr = My_find_symm_indarr(size(Vol,2));
    Zarr = My_find_symm_indarr(size(Vol,3));

    [Y, X, Z] = meshgrid(Yarr, Xarr, Zarr);

    Xshape = size(X);
    Yshape = size(Y);
    Zshape = size(Z);

    posvec = [X(:), Y(:), Z(:)]';
    rot_posvec = rotMATs*posvec;

    rotX = reshape(rot_posvec(1,:),Xshape);
    rotY = reshape(rot_posvec(2,:),Yshape);
    rotZ = reshape(rot_posvec(3,:),Zshape);

    clear posvec;
    clear rot_posvec;


    cut_size = 64;

    cutnum = floor(size(rotX,3)/cut_size);
    cutrem = mod(size(rotX,3),cut_size);

    ROTvol = zeros(size(rotX));
    
    for i=1:cutnum
        cutindar = ((i-1)*cut_size+1):i*cut_size;
        ROTtemp = interp3(Y,X,Z,Vol,rotY(:,:,cutindar),rotX(:,:,cutindar),rotZ(:,:,cutindar),'linear',0);
        ROTvol(:,:,cutindar) = ROTtemp;    
    end

    if cutrem~=0
        cutindar = (cutnum*cut_size+1):size(rotX,3);
        ROTtemp = interp3(Y,X,Z,Vol,rotY(:,:,cutindar),rotX(:,:,cutindar),rotZ(:,:,cutindar),'linear',0);
        ROTvol(:,:,cutindar) = ROTtemp;
    end

    %ROTvol = permute(ROTvol, [2 1 3]);
        
    Proj = squeeze(sum(ROTvol,3));
    
end


function dd = MatrixQuaternionRot(vector,theta)

  theta = theta*pi/180;
  vector = vector/sqrt(dot(vector,vector));
  w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
  RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
        2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
        2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];

  dd = RotM;
end
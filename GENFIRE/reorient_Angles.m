function [newAngles,convMatrix] = reorient_Angles(oldAngles,refInd,ZeroCenterFlag,varargin)

% get individual angles
Phi = oldAngles(:,1);
The = oldAngles(:,2);
Psi = oldAngles(:,3);

% get rotation matrix of refInd angle set
vector1 = [0 0 1];
rotmat1 = MatrixQuaternionRot(vector1,Phi(refInd));    

vector2 = [0 1 0];
rotmat2 = MatrixQuaternionRot(vector2,The(refInd));

vector3 = [0 0 1];
rotmat3 = MatrixQuaternionRot(vector3,Psi(refInd));

refMat = (rotmat1*rotmat2*rotmat3);

% if there is optional argument, get rotation matrix of new refInd angle
% set
if nargin > 3
    refAngle = varargin{1};
    rotmat1 = MatrixQuaternionRot(vector1,refAngle(1));      
    rotmat2 = MatrixQuaternionRot(vector2,refAngle(2));
    rotmat3 = MatrixQuaternionRot(vector3,refAngle(3));
    
    newrefMat = (rotmat1*rotmat2*rotmat3);
else
    newrefMat = eye(3);
end

% prepare empty array for output
newAngles = zeros(size(oldAngles));
convMatrix = newrefMat*refMat';

% calculate new angle set based on the reference matrix
for i=1:length(The)
  rotmat1 = MatrixQuaternionRot(vector1,Phi(i));      
  rotmat2 = MatrixQuaternionRot(vector2,The(i));
  rotmat3 = MatrixQuaternionRot(vector3,Psi(i));
  R = newrefMat*refMat'*(rotmat1*rotmat2*rotmat3);

  newAngles(i,:) = get_GENFIRE_Euler_angles_from_matrix(R,ZeroCenterFlag);
end

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
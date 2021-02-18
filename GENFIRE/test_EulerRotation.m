The = (0:5:175)' + 2;
Phi = zeros(length(The),1)+2;
Psi = zeros(length(The),1);


Ref=1;

 Phi(Ref) = Phi(Ref)-2;
%Psi(Ref) = Psi(Ref)-2;


vector1 = [0 0 1];
rotmat1 = MatrixQuaternionRot(vector1,Phi(Ref));    
% 
vector2 = [0 1 0];
rotmat2 = MatrixQuaternionRot(vector2,The(Ref));
% 
vector3 = [0 0 1];
rotmat3 = MatrixQuaternionRot(vector3,Psi(Ref));
% 

refMat = (rotmat1*rotmat2*rotmat3)';

newAngs = zeros(length(The),3);

for i=1:length(The)
rotmat1 = MatrixQuaternionRot(vector1,Phi(i));      
rotmat2 = MatrixQuaternionRot(vector2,The(i));
rotmat3 = MatrixQuaternionRot(vector3,Psi(i));

R = refMat*(rotmat1*rotmat2*rotmat3);

newAngs(i,:) = get_GENFIRE_Euler_angles_from_matrix(R);
end

[Phi newAngs(:,1) The newAngs(:,2) Psi newAngs(:,3)]
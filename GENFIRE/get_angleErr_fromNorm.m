function [TiltAngleDev, InpAngleDev] = get_angleErr_fromNorm(OriAng, NewAng)


Vec1 = [0;0;1];
Vec2 = [0;1;0];

OldVec1 = zeros(size(OriAng,1),3);
NewVec1 = zeros(size(OriAng,1),3);
OldVec2 = zeros(size(OriAng,1),3);
NewVec2 = zeros(size(OriAng,1),3);
NewVec2_proj1 = zeros(size(OriAng,1),3);

TiltAngleDev = zeros(1,size(OriAng,1));
InpAngleDev = zeros(1,size(OriAng,1));

for i=1:size(OriAng,1)

    rotmat1 = MatrixQuaternionRot([0 0 1],OriAng(i,1));
    rotmat2 = MatrixQuaternionRot([0 1 0],OriAng(i,2));   
    rotmat3 = MatrixQuaternionRot([0 0 1],OriAng(i,3));
    Matrix1 =  (rotmat1*rotmat2*rotmat3);

    % reference vectors after rotation of uncorrected angles
    OldVec1(i,:) = (Matrix1*Vec1)';
    OldVec2(i,:) = (Matrix1*Vec2)';
    
    rotmat1 = MatrixQuaternionRot([0 0 1],NewAng(i,1));
    rotmat2 = MatrixQuaternionRot([0 1 0],NewAng(i,2));   
    rotmat3 = MatrixQuaternionRot([0 0 1],NewAng(i,3));
    Matrix2 =  (rotmat1*rotmat2*rotmat3);

    % reference vectors after rotation of corrected angles
    NewVec1(i,:) = (Matrix2*Vec1)';
    NewVec2(i,:) = (Matrix2*Vec2)';
    
    % Project NewVec2 onto the plane perpendicular to OldVec1
    % so that we can calculate the in-plane rotation component
    NewVec2_OldVec1Comp = dot((OldVec1(i,:)),(NewVec2(i,:))) / (norm(OldVec1(i,:))^2) * (OldVec1(i,:));
    NewVec2_proj1(i,:) = (NewVec2(i,:))-NewVec2_OldVec1Comp;
    
    % calculate angle between normal vectors
    Dot1 = dot(OldVec1(i,:),NewVec1(i,:)) / norm(OldVec1(i,:)) / norm(NewVec1(i,:));
    if Dot1 > 1
        Dot1 = 1;
    elseif Dot1 < -1
        Dot1 = -1;
    end
    TiltAngleDev(i) = acosd( Dot1);
    
    % calculate angle between OldVec2 and projected NewVec2 onto the plane
    % normal to OldVec1. Note that OldVec2 is already peroendicular to
    % OldVec1 by definition.
    Dot2 = dot(OldVec2(i,:),NewVec2_proj1(i,:)) / norm(OldVec2(i,:)) / norm(NewVec2_proj1(i,:));
    if Dot2 > 1
        Dot2 = 1;
    elseif Dot2 < -1
        Dot2 = -1;
    end
    
    InpAngleDev(i) = acosd(Dot2 );
end
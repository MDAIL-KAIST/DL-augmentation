function [minAngle12_1, minAngle12_2] = get_CMline_bestInpAngles_LS_parallel(Pj1,Pj2,Angle1,Angle2,Range,Step)

% get rotation matrices froom the Euler angles
Matrix1 = get_Euler_matrix(Angle1);
Matrix2 = get_Euler_matrix(Angle2);

% prepare some arrays and constants
ImgSize = size(Pj1,1);
[initX, initY] = ndgrid(0,-ceil((ImgSize-1)/2):floor((ImgSize-1)/2));

% get in-plane vectors from the rotation matrices
[invnormvec1, invnormvec2] = get_commonline_vectors(Matrix1, Matrix2, [0;0;1]);

% convert vectors to angles
ComAngle1 = angle(invnormvec1(1)+invnormvec1(2)*1i) - pi/2;
ComAngle2 = angle(invnormvec2(1)+invnormvec2(2)*1i) - pi/2;

searchRangle1 = (-Range:Step:Range) + ComAngle1*180/pi;
searchRangle2 = (-Range:Step:Range) + ComAngle2*180/pi;
TotErrar = zeros(length(searchRangle1),length(searchRangle2));
parfor i=1:length(searchRangle1)
  Errar = zeros(1,length(searchRangle2))
  for j=1:length(searchRangle2)
      
    CMline1F = get_Fouriercommonline_from_angle(Pj1,searchRangle1(i), initX, initY);
    CMline2F = get_Fouriercommonline_from_angle(Pj2,searchRangle2(j), initX, initY);
    
    Errar(j) = sum(abs(CMline1F-CMline2F).^2);

  end
  TotErrar(i,:) = Errar;
end

[~,minInd] = min(TotErrar(:));
[minI,minJ] = ind2sub(size(TotErrar),minInd);

minAngle12_1 = searchRangle1(minI);
minAngle12_2 = searchRangle2(minJ);
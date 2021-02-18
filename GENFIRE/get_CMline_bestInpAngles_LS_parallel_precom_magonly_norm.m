function [minAngle12_1, minAngle12_2] = get_CMline_bestInpAngles_LS_parallel_precom_magonly_norm(PreCMline1,PreCMline2,PreAngles,Angle1,Angle2,Range)

% get rotation matrices froom the Euler angles
Matrix1 = get_Euler_matrix(Angle1);
Matrix2 = get_Euler_matrix(Angle2);

% get in-plane vectors from the rotation matrices
[invnormvec1, invnormvec2] = get_commonline_vectors(Matrix1, Matrix2, [0;0;1]);

% convert vectors to angles
ComAngle1 = angle(invnormvec1(1)+invnormvec1(2)*1i) - pi/2;
ComAngle2 = angle(invnormvec2(1)+invnormvec2(2)*1i) - pi/2;

PreAngles_moded = mod(PreAngles,360);

ComRange1_st = mod(-Range + ComAngle1*180/pi,360);
ComRange1_ed = mod(Range + ComAngle1*180/pi,360);
[~,ComRangeInd1_st] = min(abs(PreAngles_moded-ComRange1_st));
[~,ComRangeInd1_ed] = min(abs(PreAngles_moded-ComRange1_ed));

ComRange2_st = mod(-Range + ComAngle2*180/pi,360);
ComRange2_ed = mod(Range + ComAngle2*180/pi,360);
[~,ComRangeInd2_st] = min(abs(PreAngles_moded-ComRange2_st));
[~,ComRangeInd2_ed] = min(abs(PreAngles_moded-ComRange2_ed));


if ComRangeInd1_ed >= ComRangeInd1_st
  searchInds1 =  ComRangeInd1_st:ComRangeInd1_ed;
else
  searchInds1 =  [1:ComRangeInd1_ed ComRangeInd1_st:length(PreAngles_moded)];
end

if ComRangeInd2_ed >= ComRangeInd2_st
  searchInds2 =  ComRangeInd2_st:ComRangeInd2_ed;
else
  searchInds2 =  [1:ComRangeInd2_ed ComRangeInd2_st:length(PreAngles_moded)];
end

TotErrar = zeros(length(searchInds1),length(searchInds2));
parfor i=1:length(searchInds1)
  Errar = zeros(1,length(searchInds2));
  for j=1:length(searchInds2)
      
    CMline1F = PreCMline1(searchInds1(i),:);
    CMline2F = PreCMline2(searchInds2(j),:);
    
    Errar(j) = L2norm_norm(abs(CMline1F),abs(CMline2F));
  end
  TotErrar(i,:) = Errar;
end

[~,minInd] = min(TotErrar(:));
[minI,minJ] = ind2sub(size(TotErrar),minInd);

minAngle12_1 = PreAngles(searchInds1(minI));
minAngle12_2 = PreAngles(searchInds2(minJ));
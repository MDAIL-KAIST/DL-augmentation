function quatErrAr = calcQuatErrorAfterShift(p, AllAngs)

OriAng = AllAngs(:,1:3);
NewAng = AllAngs(:,4:6);




vector1 = [0 0 1];
vector2 = [0 1 0];
vector3 = [0 0 1];


  ShiftVecAng1 = p(1);
  ShiftVecAng2 = p(2);
  rotAngle = p(3);

  Mat1 = MatrixQuaternionRot([0 0 1],ShiftVecAng1);
  Mat2 = MatrixQuaternionRot([0 1 0],ShiftVecAng2);

  axisVec = Mat2*Mat1*[0 1 0]';

  finalMat = MatrixQuaternionRot(axisVec,rotAngle);

  quatErrAr = zeros(1,size(OriAng,1));

  for i=1:size(OriAng,1)
      rotmat1 = MatrixQuaternionRot(vector1,NewAng(i,1));      
      rotmat2 = MatrixQuaternionRot(vector2,NewAng(i,2));
      rotmat3 = MatrixQuaternionRot(vector3,NewAng(i,3));

      newrefMat = finalMat*(rotmat1*rotmat2*rotmat3);


      rotmat1 = MatrixQuaternionRot(vector1,OriAng(i,1));      
      rotmat2 = MatrixQuaternionRot(vector2,OriAng(i,2));
      rotmat3 = MatrixQuaternionRot(vector3,OriAng(i,3));

      OrirefMat = (rotmat1*rotmat2*rotmat3);

      quat_ori = qGetQ(OrirefMat);
      quat_new = qGetQ(newrefMat);
      quatErrAr(i) = (1 - abs(dot(quat_ori,quat_new)));
  end
    
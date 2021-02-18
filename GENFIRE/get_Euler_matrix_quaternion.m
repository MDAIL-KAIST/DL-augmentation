function Matrix = get_Euler_matrix_quaternion(Angles)
    rotmat1 = MatrixQuaternionRot([0 0 1],Angles(1));
    rotmat2 = MatrixQuaternionRot([0 1 0],Angles(2));   
    rotmat3 = MatrixQuaternionRot([0 0 1],Angles(3));
    Matrix =  (rotmat1*rotmat2*rotmat3);
end
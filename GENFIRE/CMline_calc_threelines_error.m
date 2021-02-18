function RR = CMline_calc_threelines_error(Angles1, Angles2, Angles3)%, InAngle12, InAngle13, InAngle23)

Matrix1 = get_Euler_matrix(Angles1);
Matrix2 = get_Euler_matrix(Angles2);
Matrix3 = get_Euler_matrix(Angles3);

[A_12_1, A_12_2] = get_commonline_angles(Matrix1, Matrix2, [0;0;1]);
[A_13_1, A_13_2] = get_commonline_angles(Matrix1, Matrix3, [0;0;1]);
[A_23_1, A_23_2] = get_commonline_angles(Matrix2, Matrix3, [0;0;1]);

% Err = abs(A_12_1-InAngle12(1)) + abs(A_12_2-InAngle12(2)) + ...
%       abs(A_13_1-InAngle13(1)) + abs(A_13_2-InAngle13(2)) + ...
%       abs(A_23_1-InAngle23(1)) + abs(A_23_2-InAngle23(2));

RR = mod(([A_12_1 A_12_2 A_13_1 A_13_2 A_23_1 A_23_2]+pi),pi*2)-pi;

end

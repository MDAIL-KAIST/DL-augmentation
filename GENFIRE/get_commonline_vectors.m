function [invnormvec1, invnormvec2] = get_commonline_vectors(Matrix1, Matrix2, ZeroNormal)

    % obtain rotated normal vectors
    normvec1 = Matrix1*ZeroNormal;
    normvec2 = Matrix2*ZeroNormal;
    
    % get cross correlation to get common line vector
    commonvec = cross(normvec1,normvec2);
%     commonvec = accel_cross(normvec1,normvec2);
    
    % apply inverse rotation matrix to get in-plane vectors for common line
    % for each images
    invnormvec1 = transpose(Matrix1)*commonvec;
    invnormvec2 = transpose(Matrix2)*commonvec;
end
    
    
    
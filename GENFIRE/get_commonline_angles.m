function [ComAngle1, ComAngle2] = get_commonline_angles(Matrix1, Matrix2, ZeroNormal)

    % obtain rotated normal vectors
    normvec1 = Matrix1*ZeroNormal;
    normvec2 = Matrix2*ZeroNormal;
    
    % get cross correlation to get common line vector
    commonvec = cross(normvec1,normvec2);
    
    % apply inverse rotation matrix to get in-plane vectors for common line
    % for each images
    invnormvec1 = transpose(Matrix1)*commonvec;
    invnormvec2 = transpose(Matrix2)*commonvec;
    
    % check if the vectors are indeed in-plane
    if abs(invnormvec1(3))+abs(invnormvec2(3))>0.0001
      error('bad')
    end

    % convert vectors to angles
    ComAngle1 = angle(invnormvec1(1)+invnormvec1(2)*1i) - pi/2;
    ComAngle2 = angle(invnormvec2(1)+invnormvec2(2)*1i) - pi/2;
end
    
    
    
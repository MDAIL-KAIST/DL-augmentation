% get_GENFIRE_Euler_angles_from_matrix
% Y. Yang, UCLA Physics & Astronomy
% First version date: 2016. 12. 22.

% output parameter: new_Angles: 3x1 array of phi, theta, psi angles
% input parameters: rotMAT: 3x3 rotation matrix   
%                   (optional) ZeroCenterConv: if this is true, the output
%                   theta angle will have -90 to 90 deg covention.
%                   Otherwise, 0 to 180 deg theta angle convention will be
%                   assumed,

% get_GENFIRE_Euler_angles_from_matrix calculates the Euler angles of
% phi, theta, psi in GENFIRE convention from the given rotation matrix.
% 0 to 180 degree theta angle convention will be assumed if optional
% parameter is not given. If optional parameter is true, -90 to 90 deg covention
% will be used for theta, this case assuming that phi and psi  are within
% +/- 180 degree range.




function new_Angles = get_GENFIRE_Euler_angles_from_matrix(rotMAT,varargin)

zero_tolerance = 1e-5;

if nargin > 1
    ZeroCenterConv = varargin{1};
else
    ZeroCenterConv = 0;
end

new_Angles = zeros(3,1);

if abs(rotMAT(3,3))>1
    rotMAT(3,3) = sign(rotMAT(3,3))*1;
end

new_Angles(2) = acosd(rotMAT(3,3));

if abs(sind(new_Angles(2))) > zero_tolerance
    new_Angles(1) = atan2(rotMAT(2,3),rotMAT(1,3))*180/pi;
    new_Angles(3) = atan2(rotMAT(3,2),-1*rotMAT(3,1))*180/pi;
else
    new_Angles(1) = atan2(rotMAT(3,3)*rotMAT(2,1),rotMAT(3,3)*rotMAT(1,1))*180/pi;
    new_Angles(3) = 0;
end


if (abs(new_Angles(1)) > 90 || abs(new_Angles(3)) > 90 ) && ZeroCenterConv
    MultiFactor = -1;
    
    new_Angles(2) = acosd(rotMAT(3,3))*MultiFactor;

    if abs(sind(new_Angles(2))) > zero_tolerance
        new_Angles(1) = atan2(rotMAT(2,3)*MultiFactor,rotMAT(1,3)*MultiFactor)*180/pi;
        new_Angles(3) = atan2(rotMAT(3,2)*MultiFactor,-1*rotMAT(3,1)*MultiFactor)*180/pi;
    else
        new_Angles(1) = atan2(rotMAT(2,1)*rotMAT(3,3),rotMAT(1,1)*rotMAT(3,3))*180/pi;
        new_Angles(3) = 0;
    end

end

end

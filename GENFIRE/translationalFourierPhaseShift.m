function phase_shift = translationalFourierPhaseShift(n,dx,dy,phi)
% n  - size of array (assumed square)
% dx - integer pixel shift in dim1
% dy - integer pixel shift in dim2
% phi - in plane rotation angle of line through origin

nc = round((n+1)/2); 

Y = zeros(n,1);
X = (1:n) - nc;

% just call trig functions once
c = cosd(phi);
s = sind(phi);

%calculate rotation matrix
R = [c, -s;
     s, c];

%rotate coordinates
rotkCoords = R*[X(:)';Y(:)'];

phase_shift = exp(-2j*pi*( (rotkCoords(1,:).*dx) + (rotkCoords(2,:).*dy)) / n);
end 
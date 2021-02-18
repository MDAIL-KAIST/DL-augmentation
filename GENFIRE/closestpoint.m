function x = closestpoint(n, d, p)
    % n is the vector [A,B,C] that defines the plane
    % d is the distance of the plane from the origin
    % p is the point  [P,Q,R]
    if size(p,2) == 1
        v = (d - sum(p.*n)) / sum(n.*n);
        x = p + v * n;
    else
        nr = repmat(n,[1 size(p,2)]);
        v = (d - sum(p.*nr,1)) / sum(n.*n);
        x = p + repmat(v,[3 1]) .* nr;
    end
end
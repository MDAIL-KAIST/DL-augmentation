function Err = L2norm_norm(D1,D2)
    SC = sum(D1(:).*D2(:)) / sum(D2(:).^2);
    Err = sum(abs(D1(:)-SC*D2(:)));
end
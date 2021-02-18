function Ind_shifted=My_iffshift3_ind(VolSize,Inds)

[Ind1S, Ind2S, Ind3S] = ind2sub(VolSize,Inds);
Ind1S = mod(Ind1S-1+ceil(VolSize(1)/2),VolSize(1))+1;
Ind2S = mod(Ind2S-1+ceil(VolSize(2)/2),VolSize(2))+1;
Ind3S = mod(Ind3S-1+ceil(VolSize(3)/2),VolSize(3))+1;
Ind_shifted = sub2ind(VolSize,Ind1S,Ind2S,Ind3S);
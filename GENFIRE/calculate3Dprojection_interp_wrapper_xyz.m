function calc_pj = calculate3Dprojection_interp_wrapper_xyz(obj,phi,theta,psi)
% if isempty(obj.modelK) % if empty, need to populate it. This way you only have to do this once
%     obj.modelK = my_fft(My_paddzero(obj.REFINER.refineModel,round(size(obj.REFINER.refineModel)*obj.oversampling_ratio)));
% end
% initialize arrays and get sizes
[dim1,dim2,~] = size(obj.REFINER.refineModel);
ncx           = round((dim1+1)/2);
ncy           = round((dim2+1)/2);
nc_padded1    = round((size(obj.modelK,1)+1)/2);
nc_padded2    = round((size(obj.modelK,1)+1)/2);
cropInd1      = (1:dim1)-ncx + nc_padded1;
cropInd2      = (1:dim2)-ncy + nc_padded2;

[calc_pj, kSlice] = calculate3Dprojection_interp_xyz(obj.modelK,phi,theta,psi); % get projection
% calc_pj = calc_pj(cropInd1,cropInd2).*FPmask(:,:,pj_num);
calc_pj = calc_pj(cropInd1,cropInd2);
end
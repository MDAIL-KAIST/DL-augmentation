function Vol = Func_create_voxel_gau(model, atomtype, Heights, Bfactors, AtomNumbers, volsize, Res, CropHalfWidth)
% J.Lee, KAIST, 2020

model = model ./ Res;
FHeights = Heights;
FWidths = Bfactors / pi^2 / Res^2;

if length(volsize) == 3
    x = (1:volsize(1)) - round((volsize(1)+1)/2);
    y = (1:volsize(2)) - round((volsize(2)+1)/2);
    z = (1:volsize(3)) - round((volsize(3)+1)/2);
elseif length(volsize) == 1
    x = (1:volsize(1)) - round((volsize(1)+1)/2);
    y = x;
    z = x;
else
    error('volsize should be either length 3 or length 1!')
end
    
sizeX = [length(x) length(y) length(z)];

inInd = find(model(1,:) >= min(x) & model(1,:) <= max(x) & ...
        model(2,:) >= min(y) & model(2,:) <= max(y) & ...
        model(3,:) >= min(z) & model(3,:) <= max(z));
         
calcModel = model(:,inInd);
calcAtomtype = atomtype(:,inInd);

finalvol_padded = zeros( [sizeX + (CropHalfWidth+1)*2, length(AtomNumbers)]);

cenPos = round((size(finalvol_padded)+1)/2);
cropIndRef = -CropHalfWidth:CropHalfWidth;
[cropX,cropY,cropZ] = ndgrid(cropIndRef,cropIndRef,cropIndRef);

for i=1:size(calcModel,2)
    
    currPos = calcModel(:,i) + cenPos(1:3)';
    currRndPos = round(currPos);
    
    cropInd1 = cropIndRef + currRndPos(1);
    cropInd2 = cropIndRef + currRndPos(2);
    cropInd3 = cropIndRef + currRndPos(3);
    
    CropVol = finalvol_padded(cropInd1,cropInd2,cropInd3,calcAtomtype(i));

    diffPos = currPos-currRndPos;
    
    gaussCalc = FHeights(calcAtomtype(i))*exp( -1*( (cropX-diffPos(1)).^2 + (cropY-diffPos(2)).^2 + (cropZ-diffPos(3)).^2 )/FWidths(calcAtomtype(i)) );
    
    finalvol_padded(cropInd1,cropInd2,cropInd3,calcAtomtype(i)) = CropVol + gaussCalc;
end

finalvol_summed = zeros( sizeX);


for j=1:length(AtomNumbers)
    CVol = My_stripzero(finalvol_padded(:,:,:,j),sizeX);
    finalvol_summed =finalvol_summed+CVol;
end

Vol = finalvol_summed;
    
    
    
    
    


    
    


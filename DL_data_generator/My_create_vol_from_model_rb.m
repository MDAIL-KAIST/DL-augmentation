% origin of model should be (0,0,0)
% model is 3xN array, with real unit (A, nm, etc.)
% Res is pixel size, same unit with model
% added by J.Lee (random Bfactor)

function Vol = My_create_vol_from_model_rb(model, atomtype, Heights, Bfactors, Bfactors_sigma, AtomNumbers, volsize, Res, CropHalfWidth)

model = model ./ Res;
FHeights = Heights;
FWidths = Bfactors / pi^2 / Res^2;

if ~all(Bfactors_sigma==0)
    random_bfactor=1;
end

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
    
    % random bfactor
    if random_bfactor
        endflag=1;
        while endflag
            Bfactors_tmp=Bfactors+randn([1, length(Bfactors_sigma)]).*Bfactors_sigma;

            if all(Bfactors_tmp > Bfactors/2) && all(Bfactors_tmp < Bfactors*1.5)
                endflag=0;
            end
        end
        FWidths=Bfactors_tmp / pi^2 / Res^2;
    end
    %
    
    gaussCalc = FHeights(calcAtomtype(i))*exp( -1*( (cropX-diffPos(1)).^2 + (cropY-diffPos(2)).^2 + (cropZ-diffPos(3)).^2 )/FWidths(calcAtomtype(i)) );
    
    finalvol_padded(cropInd1,cropInd2,cropInd3,calcAtomtype(i)) = CropVol + gaussCalc;
end

finalvol_summed = zeros( sizeX);

kx = 1:size(finalvol_summed,1);
ky = 1:size(finalvol_summed,2);
kz = 1:size(finalvol_summed,3);

MultF_X = 1/(length(kx)*Res);
MultF_Y = 1/(length(ky)*Res);
MultF_Z = 1/(length(kz)*Res);

CentPos = round((size(finalvol_summed)+1)/2);
[KX, KY, KZ] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y,(kz-CentPos(3))*MultF_Z);
q2 = KX.^2 + KY.^2 + KZ.^2;
clear KX KY KZ

for j=1:length(AtomNumbers)
    fa = reshape(fatom_vector( sqrt(q2),AtomNumbers(j)),sizeX);
    CVol = My_stripzero(finalvol_padded(:,:,:,j),sizeX);
    FVol = My_FFTN(CVol);
    FVol = FVol .* fa ;
    finalvol_summed =finalvol_summed+FVol;
end

Vol = real(My_IFFTN(finalvol_summed));

    
    
    
    
    
    


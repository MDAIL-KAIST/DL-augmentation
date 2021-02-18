function RandVol = My_generate_random_volume(VolSize,volumeMinMax, sphereRadiusMinMax, smoothingKernelSize)

DoneFlag = 0;
VolCenPos = round( (VolSize+1)/2);
IntLimit = VolCenPos -3; % created mask will not cover outside the IntLimit sphere

RandVol = zeros([VolSize VolSize VolSize]);

while ~DoneFlag
    SphereRad = round((rand(1)*(max(sphereRadiusMinMax) - min(sphereRadiusMinMax)))+min(sphereRadiusMinMax));
    MaxCenPos = IntLimit - SphereRad;
    SphereRadPos = (rand(1))*MaxCenPos;
    SpherePos = round( MatrixQuaternionRot([0 0 1],rand(1)*360) * MatrixQuaternionRot([0 1 0],rand(1)*360)*[SphereRadPos;0;0]);
    
    BoxInd1 = -1*SphereRad:SphereRad;
    BoxInd2 = -1*SphereRad:SphereRad;
    BoxInd3 = -1*SphereRad:SphereRad;
    BoxInd1_vol = BoxInd1 + SpherePos(1) + VolCenPos;
    BoxInd2_vol = BoxInd2 + SpherePos(2) + VolCenPos;
    BoxInd3_vol = BoxInd3 + SpherePos(3) + VolCenPos;
    
    SphereVol = zeros(SphereRad*2+1,SphereRad*2+1,SphereRad*2+1);
    [X, Y, Z] = ndgrid(BoxInd1,BoxInd2,BoxInd3);
    SphereVol(X.^2+Y.^2+Z.^2 < SphereRad^2) = 1;
    
    CurrCropVol = RandVol(BoxInd1_vol,BoxInd2_vol,BoxInd3_vol);
    CurrCropVol(SphereVol>0) = 1;
    RandVol(BoxInd1_vol,BoxInd2_vol,BoxInd3_vol) = CurrCropVol;
    
    currVol = sum(RandVol(:));
    if currVol > volumeMinMax(1)
        CC = bwconncomp(RandVol);
        if CC.NumObjects == 1 && currVol < volumeMinMax(2)
            DoneFlag = 1;
        else
            RandVol = zeros([VolSize VolSize VolSize]);
        end
    end
end
    


RandVol = smooth3(RandVol,'gaussian', smoothingKernelSize);
RandVol(RandVol>0.5) = 1;
RandVol(RandVol<1) = 0;

end
    
    
    

function CMline_F = get_Fouriercommonline_from_angle(PJ,Angle, initX, initY)   
    ImgSize = size(PJ,1);

    RotMat = [ cosd(Angle) -sind(Angle);
                sind(Angle) cosd(Angle)];

    % rotated common line coordinates
    rotVec = RotMat*[initX(:)';initY(:)'];

    % prepare arrays for DFT
    CMline_F = zeros(1,size(rotVec,2));
    
    [Xind, Yind] = ndgrid(-ceil((ImgSize-1)/2):floor((ImgSize-1)/2),-ceil((ImgSize-1)/2):floor((ImgSize-1)/2));

    N1 = ImgSize;
    N2 = ImgSize;

    % obtain Fourier common lines by DFT
    for i=1:size(rotVec,2)
        CMline_F(i) = sum(PJ(:).*exp(-1i*2*pi*(Xind(:)*rotVec(1,i)/N1 + Yind(:)*rotVec(2,i)/N2)));
    end
   
end
%%  fillInFourierGrid %%

%%inputs:
%%  projections - measured projections
%%  angles - Euler angles in the form 3xN_projections, where each projection has 3 angles in the form [phi;theta;psi]
%%  interpolationCutoffDistance - radius of sphere in which to include measured
%%  oversamplingRatio - oversampling ratio for the projections in each direction
%%      values when filling in grid. All points within this sphere will be weighted
%%      linearly by their inverse distance. 
%%  interpolationCutoffDistance - radius of interpolation kernel
%%  doCTFcorrection - flag to correct for Contrast Transfer Function (CTF) in projections, requires CTFparameters
%%  CTFparameters - structure containing defocus values and defocus angle for each projection
%%  allowMultipleGridMatches - whether or not to allow each measured datapoint to be matched to multiple grid points

%%outputs:
%%  rec - inverse FFT of the assembled Fourier grid
%%  measuredK -assembled Fourier Grid

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.




function [rec, measuredK] = fillInFourierGrid(projections, angles, particleWindowSize, oversamplingRatio, interpolationCutoffDistance, doCTFcorrection, CTFparameters, allowMultipleGridMatches)

%create empty CTF parameters if not doing CTF correction
if ~doCTFcorrection
   CTFparameters = []; 
end
if doCTFcorrection && nargin < 6
    error('GENFIRE: doCTFcorrection is turned on, but CTFparameters was not provided.\n\n')
end

%calculate padding parameters for the inputted window size
padding = round(particleWindowSize*(oversamplingRatio-1)/2);
centralPixel = size(projections,2)/2+1;
halfWindowSize = particleWindowSize/2;

%initialize array to hold measured data
kMeasured = zeros(particleWindowSize*oversamplingRatio,particleWindowSize*oversamplingRatio,size(projections,3));

tic %start clock

%get the dimension (assumed square and even) and setup the center and radius of the array size
dim1 = size(kMeasured,1);
nc = single(round((dim1+1)/2));%center pixel
n2 = single(nc-1);%radius of array

%setup the coordinates of the reciprocal slice to determine its 3D coordinates
[ky, kx] = meshgrid(-n2:n2-1,-n2:n2-1);ky = single(ky);kx = single(kx);
Q = sqrt(ky.^2+kx.^2)./n2;
kx = single(kx(:))'; ky = single(ky(:))'; %initialize coordinates of unrotate projection slice
kz = zeros(1,dim1*dim1,'single'); %0 degree rotation is a projection onto the X-Y plane, so all points have kz=0;

%check for the presence of some of the CTF correction options and set defaults if they are absent
if doCTFcorrection
    if isfield(CTFparameters,'CTFThrowOutThreshhold')
        CTFThrowOutThreshhold = CTFparameters(1).CTFThrowOutThreshhold; %value below which to not grid points that were suppressed by the CTF
    else
        CTFThrowOutThreshhold = 0.05;%default value
    end
    if isfield(CTFparameters,'ignore_first_peak')
       ignore_first_peak =  CTFparameters(1).ignore_first_peak;
    else
        ignore_first_peak = 0;
    end  

for projNum = 1:size(projections,3);
    %get Contrast Transfer Function (CTF)
    pjK = projections(:,:,projNum);
    centralPixelK = size(pjK,2)/2+1;
    
    %crop out the appropriate window
    pjK = pjK(centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1,centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1);%window projection

    pjK = my_fft(padarray(pjK,[padding padding 0]));%pad and take FFT
    
    %get the CTF
    [CTF, gamma] = ctf_correction(pjK,CTFparameters(projNum).defocusU,CTFparameters(projNum).defocusV,CTFparameters(projNum).defocusAngle,ignore_first_peak);%get CTF
    if CTFparameters(projNum).phaseFlip %this should always be on unless your projections have already been CTF corrected elsewhere
        pjK(CTF<0) = -1*pjK(CTF<0);%phase flip
    end
    
    if CTFparameters(projNum).correctAmplitudesWithWienerFilter
    	
    	%get dimensions of the CTF array
        dim1_2 = size(CTF,1);
        nc2 = single(round((dim1_2+1)/2));%center pixel
        n22 = single(nc2-1);%radius of array

		%reciprocal indices
        [ky2, kx2] = meshgrid(-n22:n22-1,-n22:n22-1);ky2 = single(ky2);kx2 = single(kx2);
        Q2 = sqrt(ky2.^2+kx2.^2)./n22;
        
        SSNR = ones(size(Q2));%initialize SSNR map
        %interpolate the SSNR array from the provided values of the SSNR per frequency shell
        SSNR(:) = interp1(linspace(0,1+1e-10,size(CTFparameters(projNum).SSNR,2)),CTFparameters(projNum).SSNR,Q2(:),'linear');%make weighting map from average FRC
        SSNR(isnan(SSNR)) = 0;
        wienerFilter = abs(CTF)./(abs(CTF).^2+(1./SSNR));%construct Wiener filter for CTF amplitude correction
        pjK = pjK.*wienerFilter; 
    elseif CTFparameters(projNum).multiplyByCTFabs%multiplying by CTF boosts SNR and is most useful for datasets that are extremely noisy
        pjK = pjK.*abs(CTF); 
    end
    
    if CTFThrowOutThreshhold>0 %recalculate CTF at new array size for throwing out values that were near CTF 0 crossover
        pjK(abs(CTF)<CTFThrowOutThreshhold & (gamma>(pi/2))) = -999;%flag values where CTF was near 0 to ignore for gridding, but ignore out to first peak
    end
    
    kMeasured(:,:,projNum) = pjK;   
end
    
    if CTFThrowOutThreshhold > 0     %flag values below the where the CTF was smaller than the CTFThrowOutThreshhold
        for projNum = 1:size(projections,3);
            pjK = projections(centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1,centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1,projNum);
            pjK = my_fft(padarray(pjK,[padding padding 0]));
            CTF = ctf_correction(pjK,CTFparameters(projNum).defocusU,CTFparameters(projNum).defocusV,CTFparameters(projNum).defocusAngle,ignore_first_peak);%get CTF
            pjK(abs(CTF)<CTFThrowOutThreshhold) = -999;%flag values where CTF was near 0 to ignore for gridding
            kMeasured(:,:,projNum) = pjK;
        end  
    end
else
    %otherwise, add the projection to the stack of data with no further corrections
    for projNum = 1:size(projections,3);
        kMeasured(:,:,projNum) = my_fft(padarray(projections(centralPixel-halfWindowSize:centralPixel+halfWindowSize-1,centralPixel-halfWindowSize:centralPixel+halfWindowSize-1,projNum),[padding padding  0]));
    end  
end

clear projections

%initialize arrays to contain coordinates
measuredX = zeros(1,size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'single');
measuredY = zeros(1,size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'single');
measuredZ = zeros(1,size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'single');


for projNum = 1:size(kMeasured,3);
phi = angles(projNum,1);
theta = angles(projNum,2);
psi = angles(projNum,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  GENFIRE/RELION/XMIPP/FREALIGN/EMAN Euler angle convention:
% % 
R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
      -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
      sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];

rotkCoords = R'*[kx;ky;kz];%rotate coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

measuredX(:,:,projNum) = rotkCoords(1,:);%rotated X
measuredY(:,:,projNum) = rotkCoords(2,:);%rotated Y
measuredZ(:,:,projNum) = rotkCoords(3,:);%rotated Z
end

%reshape to simplify
measuredX = reshape(measuredX,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
measuredY = reshape(measuredY,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
measuredZ = reshape(measuredZ,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
kMeasured = reshape(kMeasured,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
badInd = find(kMeasured==-999);%delete values that are flagged as bad
measuredX(badInd) = [];
measuredY(badInd) = [];
measuredZ(badInd) = [];
kMeasured(badInd) = [];

masterInd = [];%masterInd will be a large list of the grid indices
masterVals = [];%complex values to include in weighted averaging for those grid points
masterDistances = [];%distance from measured value to grid point
% masterConfidenceWeights = [];

if allowMultipleGridMatches
    shiftMax = round(interpolationCutoffDistance);
else
    shiftMax = 0;
end
%The nearest grid point to a measured value can be found by rounding, but
%there can be more than one grid point within the cutoff sphere, so must
%search locally for other possibilities. However in practice I have found 
%this search can slow the program down greatly, without significant change
%in final result. Even searching 1 voxel in either direction increases the
%number of calculations by 3^3 = 27; For this reason I have set shiftMax = 0 and
%just assign values to their closest voxel.

for Yshift = -shiftMax:shiftMax 
   for Xshift = -shiftMax:shiftMax
       for Zshift = -shiftMax:shiftMax
            tmpX = (round(measuredX)+Xshift); % apply shift
            tmpY = (round(measuredY)+Yshift);
            tmpZ = (round(measuredZ)+Zshift);
            tmpVals = kMeasured;
            distances = sqrt(abs(measuredX-tmpX).^2+abs(measuredY-tmpY).^2+abs(measuredZ-tmpZ).^2); %compute distance to nearest voxel
            tmpY = tmpY+nc; %shift origin
            tmpZ = tmpZ+nc;
            tmpX = tmpX+nc;
            goodInd = (~(tmpX>dim1|tmpX<1|tmpY>dim1|tmpY<1|tmpZ>dim1|tmpZ<1)) & distances<=interpolationCutoffDistance;%find candidate values
            masterInd = [masterInd sub2ind([dim1 dim1 dim1],tmpX(goodInd),tmpY(goodInd),tmpZ(goodInd))]; %append values to lists
            masterVals = [masterVals tmpVals(goodInd)];
            masterDistances = [masterDistances distances(goodInd)];

       end
   end
end
   
clear measuredX
clear measuredY
clear measuredZ
clear confidenceWeights

% Now that we have a list of the complex values to grid, their coordinates, 
% and their distances from the nearest voxel, we want to reorganize the
% data so that all values matched to a given voxel are in the same place,
% so that the weighted sum can be computed. The number of values matched to
% each voxel can vary, and although one could use cell arrays for this
% purpose, they are quite slow. Instead, one can simply sort the indices,
% and then find the unique values by looking at the difference in
% consecutive elements. 

masterDistances = masterDistances + 1e-5;
masterDistances(masterDistances>0) = 1 ./ masterDistances(masterDistances>0);
masterDistances(isnan(masterDistances)) = 0;

measuredK = accumarray(masterInd',masterVals.*masterDistances,[dim1^3 1]);
sumWeights = accumarray(masterInd',masterDistances,[dim1^3 1]);
measuredK(sumWeights>0) = measuredK(sumWeights>0) ./ sumWeights(sumWeights>0);
measuredK = reshape(measuredK,[dim1 dim1 dim1]);
measuredK = hermitianSymmetrize(measuredK);

rec = real(my_ifft(measuredK));
timeTakenToFillInGrid = toc;
timeTakenToFillInGrid = round(10*timeTakenToFillInGrid)./10;
fprintf('GENFIRE: Fourier grid assembled in %.12g seconds.\n\n',timeTakenToFillInGrid);

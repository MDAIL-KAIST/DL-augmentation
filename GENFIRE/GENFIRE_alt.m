%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                         %%
%%                        Welcome to GENFIRE!                              %%
%%           GENeralized Fourier Iterative REconstruction                  %%
%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Author: Alan (AJ) Pryor, Jr.
%% email:  apryor6@gmail.com
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GENFIRE_alt(filename_Projections,filename_Angles,filename_Support,filename_Results)
addpath ./src/
addpath ./data/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLOT_YN = 1;

GENFIRE = GENFIRE_Reconstructor();

%%% See the README for description of parameters

GENFIRE.filename_Projections = filename_Projections;
GENFIRE.filename_Angles = filename_Angles;
GENFIRE.filename_Support = filename_Support;
% filename_InitialModel = '';
GENFIRE.filename_Results = filename_Results;

GENFIRE.numIterations = 50; 
GENFIRE.pixelSize = .5; 
GENFIRE.oversamplingRatio =3;
GENFIRE.griddingMethod = 2; 
GENFIRE.allowMultipleGridMatches = 1;
GENFIRE.constraintEnforcementMode = 1; 
GENFIRE.interpolationCutoffDistance =.7; 
GENFIRE.constraintPositivity = 1;
GENFIRE.constraintSupport = 1;
GENFIRE.ComputeFourierShellCorrelation = 0; 
GENFIRE.numBins = 50;
GENFIRE.percentValuesForRfree = 0.05;
GENFIRE.numBinsRfree = 35;
GENFIRE.doCTFcorrection = 0;
GENFIRE.CTFThrowOutThreshhold = 0;
%GENFIRE.particleWindowSize = [];
%GENFIRE.phaseErrorSigmaTolerance = [];
%GENFIRE.vector1 = [0 0 1];
%GENFIRE.vector2 = [0 1 0];
%GENFIRE.vector3 = [0 0 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin GENFIRE

% Read data from files. The data can also be directly be fed into GENFIRE
% Class, not necessarily from file
GENFIRE = readFiles(GENFIRE);

% Check if the data is consistent, and prepare reconstruction variables
% based on given parameters
GENFIRE = CheckPrepareData(GENFIRE);

% Run GENFIRE gridding
GENFIRE = runGridding(GENFIRE); 

% Run FSC if flag is on
if GENFIRE.ComputeFourierShellCorrelation
  
    % Do not calculate R_free for FSC calculation
    calculate_Rfree_ori = GENFIRE.calculate_Rfree;
    GENFIRE.calculate_Rfree= 0;
    
    GENFIRE = runFSC(GENFIRE);
    
    % set the R_free flag back go original
    GENFIRE.calculate_Rfree=calculate_Rfree_ori;
  
    if PLOT_YN == 1
        figure, plot(GENFIRE.spatialFrequency,GENFIRE.FSC,'k','LineWidth',3)
        set(gcf,'color','white')
        title('FSC between independent half reconstructions','FontSize',16)
        xlabel('Spatial Frequency','FontSize',14)
        ylabel('Correlation Coefficient','FontSize',14)
    
    end
end
    
%run reconstruction

GENFIRE = GENFIRE.reconstruct();

if PLOT_YN == 1
  ncX = round((GENFIRE.Dim1+1)/2);
  ncY = round((GENFIRE.Dim2+1)/2);
  ncZ = round((GENFIRE.Dim1+1)/2);

  subplot(2,3,4), imagesc(squeeze(sum(GENFIRE.reconstruction,1))),title('GENFIRE projection 1'), axis image
  subplot(2,3,5), imagesc(squeeze(sum(GENFIRE.reconstruction,2))),title('GENFIRE projection 2'), axis image
  subplot(2,3,6), imagesc(squeeze(sum(GENFIRE.reconstruction,3))),title('GENFIRE projection 3'), axis image
  subplot(2,3,1), imagesc(squeeze(sum(GENFIRE.recIFFT,1))),title('before iteration projection 1'), axis image
  subplot(2,3,2), imagesc(squeeze(sum(GENFIRE.recIFFT,2))),title('before iteration projection 2'), axis image
  subplot(2,3,3), imagesc(squeeze(sum(GENFIRE.recIFFT,3))),title('before iteration projection 3'), axis image


  figure,
  subplot(2,3,4), imagesc(squeeze(GENFIRE.reconstruction(ncX,:,:))),title('GENFIRE slice 1'), axis image
  subplot(2,3,5), imagesc(squeeze(GENFIRE.reconstruction(:,ncY,:))),title('GENFIRE slice 2'), axis image
  subplot(2,3,6), imagesc(squeeze(GENFIRE.reconstruction(:,:,ncZ))),title('GENFIRE slice 3'), axis image
  subplot(2,3,1), imagesc(squeeze(GENFIRE.recIFFT(ncX,:,:))),title('before iteration slice 1'), axis image
  subplot(2,3,2), imagesc(squeeze(GENFIRE.recIFFT(:,ncY,:))),title('before iteration slice 2'), axis image
  subplot(2,3,3), imagesc(squeeze(GENFIRE.recIFFT(:,:,ncZ))),title('before iteration slice 3'), axis image
end


GENFIRE = ClearCalcVariables(GENFIRE);
    
SaveResults(GENFIRE);
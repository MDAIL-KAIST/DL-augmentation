
% This script calculates projections and reconstructs with perturbed
% angles which are then recovered by refinement

%addpath ./src/
model_filename   = '../models/FePt_Vol_64.mat';
pj_filename      = '../data/FprojectionsT.mat';
angle_filename   = '../data/FanglesT.mat';
support_filename = '../data/FsupportT.mat';

random_seed     = 42; % for reproducibility
tomo_angles     = -90:5:85; % true tilt angles
OR              = 4; % oversampling ratio used for padding model in 3D for projection calculation
jitter          = 4; % euler angles will be perturbed by up to += jitter/2
num_refinements = 1; % number of refinement loops

tomo_full_angles = [zeros(length(tomo_angles),1) tomo_angles' zeros(length(tomo_angles),1)];
%% calculate projections
model = importdata(model_filename);
[dim,~,~] = size(model);
nc = round((dim+1)/2); n2 = nc-1;
padding = (OR-1)*dim;
modelK = my_fft(padarray(model,[padding, padding, padding]));
nc_padded = round((size(modelK,1)+1)/2);
num_proj = length(tomo_angles);
pj = zeros(dim,dim,num_proj);


%%
exp_angles_jittered = [zeros(length(tomo_angles),1) tomo_angles' zeros(length(tomo_angles),1)];
jittered_ind = 1:36;
rng(42);
exp_angles_jittered(jittered_ind,1) = exp_angles_jittered(jittered_ind,1) + jitter*rand(length(jittered_ind),1) - jitter/2;
rng(43);
exp_angles_jittered(jittered_ind,2) = exp_angles_jittered(jittered_ind,2) + jitter*rand(length(jittered_ind),1) - jitter/2;
rng(44);
exp_angles_jittered(jittered_ind,3) = exp_angles_jittered(jittered_ind,3) + jitter*rand(length(jittered_ind),1) - jitter/2;


%%
parfor pj_num = 1:num_proj
    if mod(pj_num,1)==0
        fprintf('Calculating projection #%d/%d\n',pj_num,num_proj)
    end
   tmp = calculate3Dprojection_interp(modelK,exp_angles_jittered(pj_num,1),exp_angles_jittered(pj_num,2),exp_angles_jittered(pj_num,3)); 
   pj(:,:,pj_num) = tmp(nc_padded-n2:nc_padded+n2-1,nc_padded-n2:nc_padded+n2-1);
end

%% perturb angles

support = ones(dim,dim,dim);

%%
save(pj_filename,'pj')
save(angle_filename,'exp_angles_jittered')
save(support_filename,'support')

%%
Ref = 18;

Phi = exp_angles_jittered(:,1);
The = exp_angles_jittered(:,2);
Psi = exp_angles_jittered(:,3);

vector1 = [0 0 1];
rotmat1 = MatrixQuaternionRot(vector1,Phi(Ref));    
% 
vector2 = [0 1 0];
rotmat2 = MatrixQuaternionRot(vector2,The(Ref));
% 
vector3 = [0 0 1];
rotmat3 = MatrixQuaternionRot(vector3,Psi(Ref));
% 

refMat = (rotmat1*rotmat2*rotmat3)';

newAngs = zeros(length(The),3);

for i=1:length(The)
rotmat1 = MatrixQuaternionRot(vector1,Phi(i));      
rotmat2 = MatrixQuaternionRot(vector2,The(i));
rotmat3 = MatrixQuaternionRot(vector3,Psi(i));

R = refMat*(rotmat1*rotmat2*rotmat3);

newAngs(i,:) = get_GENFIRE_Euler_angles_from_matrix(R);
end




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


%addpath ./src/
addpath ./data/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLOT_YN = 1;

GENFIRE = GENFIRE_Class();

%%% See the README for description of parameters

GENFIRE.filename_Projections = '../data/FprojectionsT.mat';
GENFIRE.filename_Angles = '../data/FanglesT.mat';
%GENFIRE.filename_Support = './data/Emptysupport.mat'; 
% filename_InitialModel = '';
GENFIRE.filename_Results = './results/GENFIRE_recEulerTest.mat';

GENFIRE.numIterations = 50; 
GENFIRE.pixelSize = .5; 
GENFIRE.oversamplingRatio =3;
GENFIRE.griddingMethod = 2; 
GENFIRE.allowMultipleGridMatches = 1;
GENFIRE.constraintEnforcementMode = 3; 
GENFIRE.interpolationCutoffDistance =.3; 
GENFIRE.constraintPositivity = 1;
GENFIRE.constraintSupport = 1;
GENFIRE.ComputeFourierShellCorrelation = 0; 
GENFIRE.numBins = 50;
GENFIRE.percentValuesForRfree = 0.0;
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

GENFIRE.vector2=[0 1 0];
% Check if the data is consistent, and prepare reconstruction variables
% based on given parameters
%GENFIRE=GENFIRE.set_parameters('InputAngles',exp_angles_jittered);
GENFIRE=GENFIRE.set_parameters('InputAngles',newAngs);


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

%%
%run reconstruction
GENFIRE=GENFIRE.set_parameters('calculate_Rfree',0);
GENFIRE=GENFIRE.set_parameters('RandomSeed',50);

tic
GENFIRE2 = runGENFIREiteration_GENFIREClass(GENFIRE);
toc
% tic
% GENFIRE2 = runGENFIREiteration_GENFIREClass_old(GENFIRE);
% toc


%%
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
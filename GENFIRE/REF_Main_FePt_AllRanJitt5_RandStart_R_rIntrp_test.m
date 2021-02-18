% This script calculates projections and reconstructs with perturbed
% angles which are then recovered by refinement

%addpath ./src/
model_filename   = '../models/FePt_Vol_64.mat';
pj_filename      = '../data/RealPJ_AllRanJitt5_RI_test.mat';
angle_filename   = '../data/Angles_AllRanJitt5_RI_test.mat';
support_filename = '../data/Fsupport_RI_test.mat';

random_seed     = 101; % for reproducibility
tomo_angles     = -81:9:81; % true tilt angles
OR              = 4; % oversampling ratio used for padding model in 3D for projection calculation
jitter          = 8; % euler angles will be perturbed by up to += jitter/2
jitter_tomo          = 2; % euler angles will be perturbed by up to += jitter/2

num_refinements = 7; % number of refinement loops

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
jittered_ind = 1:length(tomo_angles);
rng(42);
exp_angles_jittered(jittered_ind,1) = exp_angles_jittered(jittered_ind,1) + jitter*rand(length(jittered_ind),1) - jitter/2;
rng(443);
exp_angles_jittered(jittered_ind,2) = exp_angles_jittered(jittered_ind,2) + jitter*rand(length(jittered_ind),1) - jitter/2 ;
rng(244);
exp_angles_jittered(jittered_ind,3) = exp_angles_jittered(jittered_ind,3) + jitter*rand(length(jittered_ind),1) - jitter/2;


%%
parfor pj_num = 1:num_proj
    if mod(pj_num,1)==0
        fprintf('Calculating projection #%d/%d\n',pj_num,num_proj)
    end
%    tmp = calculate3Dprojection_interp(modelK,exp_angles_jittered(pj_num,1),exp_angles_jittered(pj_num,2),exp_angles_jittered(pj_num,3)); 
%    pj(:,:,pj_num) = tmp(nc_padded-n2:nc_padded+n2-1,nc_padded-n2:nc_padded+n2-1);
    tmp = calculate3Dprojection(model,exp_angles_jittered(pj_num,1),exp_angles_jittered(pj_num,2),exp_angles_jittered(pj_num,3)); 
    pj(:,:,pj_num) = tmp;%(nc_padded-n2:nc_padded+n2-1,nc_padded-n2:nc_padded+n2-1);
end

%% perturb angles

support = ones(dim,dim,dim);

%%
save(pj_filename,'pj')
tomo_full_angles(:,3) = tomo_full_angles(:,3) +0;
save(angle_filename,'tomo_full_angles')
save(support_filename,'support')

%%
NameCell = {'Zero','Puturb1','Puturb2','Puturb3'};

for Loopnum = 3:length(NameCell);

resultFileName = sprintf('./results/RefineResult_AllRanJutt5_RandStart_R_RealInterp_%s_test.mat',NameCell{Loopnum});

%% run refinement

% initilize REFINEMENT object with RECONSTRUCTOR, in this case with GENFIRE
REFINEMENT = REFINE_Class('GENFIRE');

% set REFINEMENT parameters
REFINEMENT = REFINEMENT.set_parameters('ang_search_range',-5:0.25:5);
REFINEMENT = REFINEMENT.set_parameters('ang_search_range_cell',{-5:1:5,-3:1:3,-2:1:2,-1:0.25:1,-1:0.1:1,-1:0.1:1,-1:0.1:1});
REFINEMENT = REFINEMENT.set_parameters('compare_func',@alignByRFactor_M_subpixel);
REFINEMENT = REFINEMENT.set_parameters('forwardProjection_func',@calculate3Dprojection_RealSpaceinterp,'RealProjection',1);
% %REFINEMENT = REFINEMENT.set_parameters('compare_func',@alignByNormXCorr);
REFINEMENT = REFINEMENT.set_parameters('Rmethod',1);
REFINEMENT = REFINEMENT.set_parameters('maximize',false);
REFINEMENT = REFINEMENT.set_parameters('use_parallel',true);

REFINEMENT = REFINEMENT.set_parameters('RefineReferenceAngleInd',10);
REFINEMENT = REFINEMENT.set_parameters('RefineReferenceAngletoSet',[0 0 0]);
REFINEMENT = REFINEMENT.set_parameters('RefineZeroCenterFlag',1);

REFINEMENT = REFINEMENT.set_parameters('FullEvolutionRecord',1);
REFINEMENT = REFINEMENT.set_parameters('oversampling_ratio',3);

REFINEMENT = REFINEMENT.set_parameters('num_refinements',num_refinements);

% if Loopnum == 2
%  REFINEMENT = REFINEMENT.set_parameters('forwardProjection_func',@calculate3Dprojection_RealSpaceinterp,'RealProjection',1);
% else
%   if Loopnum == 1
%     REFINEMENT = REFINEMENT.set_parameters('oversampling_ratio',2);
%   elseif Loopnum == 3
%     REFINEMENT = REFINEMENT.set_parameters('oversampling_ratio',4);
%   elseif Loopnum == 4
%     REFINEMENT = REFINEMENT.set_parameters('oversampling_ratio',6);
%   end
% end



% set RECONSTRUCTOR parameters (in this case GENFIRE)
REFINEMENT.RECONSTRUCTOR.filename_Projections = pj_filename;
REFINEMENT.RECONSTRUCTOR.filename_Angles = angle_filename;
REFINEMENT.RECONSTRUCTOR.filename_Support = support_filename; 
%REFINEMENT.RECONSTRUCTORfilename_InitialModel = '';
REFINEMENT.RECONSTRUCTOR.filename_Results = './results/REFINEMENT_finalREconstruction.mat';

REFINEMENT.RECONSTRUCTOR.numIterations = 50; 
REFINEMENT.RECONSTRUCTOR.pixelSize = .5; 
REFINEMENT.RECONSTRUCTOR.oversamplingRatio =3;
REFINEMENT.RECONSTRUCTOR.griddingMethod = 2; 
REFINEMENT.RECONSTRUCTOR.allowMultipleGridMatches = 1;
REFINEMENT.RECONSTRUCTOR.constraintEnforcementMode = 3; 
REFINEMENT.RECONSTRUCTOR.interpolationCutoffDistance =.3; 
REFINEMENT.RECONSTRUCTOR.constraintPositivity = 1;
REFINEMENT.RECONSTRUCTOR.constraintSupport = 1;
REFINEMENT.RECONSTRUCTOR.ComputeFourierShellCorrelation = 0; 
REFINEMENT.RECONSTRUCTOR.numBins = 50;
REFINEMENT.RECONSTRUCTOR.percentValuesForRfree = 0.00;
REFINEMENT.RECONSTRUCTOR.numBinsRfree = 35;
REFINEMENT.RECONSTRUCTOR.doCTFcorrection = 0;
REFINEMENT.RECONSTRUCTOR.CTFThrowOutThreshhold = 0;



% initilize RECONSTRUCTOR, read data, and check data validity
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.readFiles();
if Loopnum > 1
%   rng('shuffle')
%   tomo_full_angles_jittered = tomo_full_angles + jitter_tomo*rand(size(tomo_full_angles)) - jitter_tomo/2;
%   REFINEMENT.RECONSTRUCTOR.InputAngles = tomo_full_angles_jittered;
resultFileName = sprintf('./results/RefineResult_AllRanJutt5_RandStart_R_RealInterp_%s.mat',NameCell{Loopnum});

RR = importdata(resultFileName);

 AngleEvolution = RR.REFINEMENT.AngleEvolution;
 REFINEMENT.RECONSTRUCTOR.InputAngles =AngleEvolution(:,:,1);
end
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.CheckPrepareData();


% import projections and angles from RECONSTRUCTOR to REFINEMENT
REFINEMENT = REFINEMENT.get_projections_from_RECONSTRUCTOR();  
REFINEMENT = REFINEMENT.get_angles_from_RECONSTRUCTOR();


% run refinement
tic

REFINEMENT = REFINEMENT.refineControl_REFINEClass();

toc


% get refined result
refined_angles = REFINEMENT.refineAngles;

true_angles = exp_angles_jittered;%[zeros(length(tomo_angles),1), tomo_angles', zeros(length(tomo_angles),1)];


REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.ClearCalcVariables();

save(resultFileName,'refined_angles','true_angles','REFINEMENT');

end

%%
OO = importdata('tempIter2.mat');
OO = OO.updateProjections_REFINEClass();

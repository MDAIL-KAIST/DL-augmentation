% This script calculates projections and reconstructs with perturbed
% angles which are then recovered by refinement

%addpath ./src/
model_filename   = '../models/vesicle.mat';
pj_filename      = '../data/projections.mat';
angle_filename   = '../data/angles.mat';
support_filename = '../data/support.mat';

random_seed     = 42; % for reproducibility
tomo_angles     = 0:5:180; % true tilt angles
OR              = 4; % oversampling ratio used for padding model in 3D for projection calculation
jitter          = 4; % euler angles will be perturbed by up to += jitter/2
num_refinements = 1; % number of refinement loops

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
parfor pj_num = 1:num_proj
    if mod(pj_num,1)==0
        fprintf('Calculating projection #%d/%d\n',pj_num,num_proj)
    end
   tmp = calculate3Dprojection_interp(modelK,0,tomo_angles(pj_num),0); 
   pj(:,:,pj_num) = tmp(nc_padded-n2:nc_padded+n2-1,nc_padded-n2:nc_padded+n2-1);
end

%% perturb angles
euler_angles = zeros(num_proj,3);
rng(random_seed);
euler_angles(:,2) = tomo_angles + round(jitter*rand(length(tomo_angles),1)' - jitter/2);
support = ones(dim,dim,dim);

%%
save(pj_filename,'pj')
save(angle_filename,'euler_angles')
save(support_filename,'support')


%% run refinement

% initilize REFINEMENT object with RECONSTRUCTOR, in this case with GENFIRE
REFINEMENT = REFINE_Class('GENFIRE');

% set REFINEMENT parameters
REFINEMENT = REFINEMENT.set_parameters('ang_search_range',-2:2,'compare_func',@alignByNormXCorr, 'maximize',true,'use_parallel',true);
REFINEMENT = REFINEMENT.set_parameters('FullEvolutionRecord',1);
%REFINEMENT = REFINEMENT.set_parameters('forwardProjection_func',@calculate3Dprojection_RealSpaceinterp,'RealProjection',1);
REFINEMENT = REFINEMENT.set_parameters('oversampling_ratio',4);



% set RECONSTRUCTOR parameters (in this case GENFIRE)
REFINEMENT.RECONSTRUCTOR.filename_Projections = pj_filename;
REFINEMENT.RECONSTRUCTOR.filename_Angles = angle_filename;
REFINEMENT.RECONSTRUCTOR.filename_Support = support_filename; 
%REFINEMENT.RECONSTRUCTORfilename_InitialModel = '';
REFINEMENT.RECONSTRUCTOR.filename_Results = './results/REFINEMENT_finalREconstruction.mat';

REFINEMENT.RECONSTRUCTOR.numIterations = 50; 
REFINEMENT.RECONSTRUCTOR.pixelSize = .5; 
REFINEMENT.RECONSTRUCTOR.oversamplingRatio =4;
REFINEMENT.RECONSTRUCTOR.griddingMethod = 1; 
REFINEMENT.RECONSTRUCTOR.allowMultipleGridMatches = 1;
REFINEMENT.RECONSTRUCTOR.constraintEnforcementMode = 3; 
REFINEMENT.RECONSTRUCTOR.interpolationCutoffDistance =.7; 
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
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.CheckPrepareData();


% import projections and angles from RECONSTRUCTOR to REFINEMENT
REFINEMENT = REFINEMENT.get_projections_from_RECONSTRUCTOR();  
REFINEMENT = REFINEMENT.get_angles_from_RECONSTRUCTOR();


% run refinement
REFINEMENT = REFINEMENT.refineControl_REFINEClass(num_refinements);


% get refined result
refined_angles = REFINEMENT.refineAngles;

true_angles = [zeros(length(tomo_angles),1), tomo_angles', zeros(length(tomo_angles),1)];

% compare results
plotResults(true_angles,refined_angles);
figure;plot(true_angles(:,2)-refined_angles(:,2));
figure;plot(true_angles(:,2)-euler_angles(:,2));


REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.ClearCalcVariables();

save('REFINERESULT_new_OR4.mat','refined_angles','true_angles','REFINEMENT');


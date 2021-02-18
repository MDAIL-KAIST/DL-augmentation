clear all
clc
%%
% Author: J.Lee, KAIST (Korea), 2021.
% Y.Yang, Multi-Dimensional Atomic Imaging Lab, KAIST
% Generating input/target data for deep learning augmentation.
% Amorphous structure version

%% parameter setting
DLDgen = DL_Dataset_Generator_Class();

% common parameters
DLDgen.Res = 0.3570; % pixel size
DLDgen.volsize = 48; % volume size
DLDgen.Number_data_start = 1; % start number
DLDgen.Number_data_final = 1; % end number
DLDgen.Atomnumbers = [78]; % atomnumber (ex, 78 -> Pt, 26 -> Fe)

% input & target saving location
DLDgen.TARGET_prefix = 'Pt_target'; 
DLDgen.TARGET_PATH = './Pt_target';
DLDgen.INPUT_prefix = 'Pt_input';
DLDgen.INPUT_PATH = './Pt_input';

%%% target data %%%
% target data parameters
DLDgen.target_gen_YN = 1; % create target data on/off
% 
% random parameters for target data
DLDgen.target_Random_shape = 1; % crop shape randomly on/off
%
DLDgen.target_atomdensity = 4/3.912^3 ; % target atomdensity
DLDgen.target_minimum_distance = 2.0 ; % minimum distance between atoms (Angstron)
DLDgen.target_getting_percent = 0.99; % generating atoms until getting target percent of atomdensity
%
DLDgen.target_gen_Gaus_YN = 1; % create target data (ground truth) on/off  
DLDgen.target_Height = 78; % intensity scale for target data
DLDgen.target_Bfactors = 4; % Bfactor for input data
DLDgen.target_Atomnumbers = DLDgen.Atomnumbers;
DLDgen.target_CropHalfWidth = 4; % crop boxsize for generating gaussian function

%%% input data %%%
% input data parameters
DLDgen.input_data_YN = 1; % create input data on/off
%
DLDgen.input_Random_missing_wedge = 0; % adjust missing wedge random angle
DLDgen.input_theta_angle_range = -65:65/10:65; % set tilt angle range
DLDgen.input_Angle_noise_option = 2; % 0: no, 1: theta only, 2: all
DLDgen.input_Angle_noise_range = 0.6; % angle error (0.6 means that -0.3 to 0.3 degree)

DLDgen.input_Height = [10]; % intensity scale for input data
DLDgen.input_Bfactors = [5]; % mean Bfactor for input data
DLDgen.input_Bfactors_sigma = [1]; % std Bfactor for input data
DLDgen.input_Atomnumbers = DLDgen.Atomnumbers;
DLDgen.input_CropHalfWidth = 8; % crop boxsize for generating gaussian function
DLDgen.input_Poisson_noise = 1; % Poisson noise on/off     

% GENFIRE parameter (Reconstruction)
DLDgen.GENFIRE.numIterations = 100; % number of iterations
DLDgen.GENFIRE.oversamplingRatio = 2; 
DLDgen.GENFIRE.griddingMethod = 1;  % 1 for FFT, 2 for DFT
DLDgen.GENFIRE.allowMultipleGridMatches = 1; % only for FFT, this parameter will be ignored for DFT interpolation
DLDgen.GENFIRE.constraintEnforcementMode = 3; % Determines the Fourier expansion / supression behavior, "3" turns it off.
DLDgen.GENFIRE.interpolationCutoffDistance = 0.3; % interpolation radius 
DLDgen.GENFIRE.constraintPositivity = 1; 
DLDgen.GENFIRE.constraintSupport = 1;
%DLDgen.GENFIRE.doCTFcorrection = 0;
%DLDgen.GENFIRE.CTFThrowOutThreshhold = 0;
%DLDgen.GENFIRE.calculate_Rfree = 0; % flag for calculation of Rfree
%DLDgen.GENFIRE.DFT_doGPU = 0; % flag for using GPU for gridding
%DLDgen.GENFIRE.vector1 = [0 0 1];
%DLDgen.GENFIRE.vector2 = [0 1 0];
DLDgen.GENFIRE.vector3 = [1 0 0];
%%
%% running code
DLDgen.MAKEDIR(); % generate folder
DLDgen = DLDgen.GEN_target_data_amorphous(); % generate target data (Ground truth) (run first)
DLDgen = DLDgen.GEN_input_data(); % generate input data (based on Ground truth)


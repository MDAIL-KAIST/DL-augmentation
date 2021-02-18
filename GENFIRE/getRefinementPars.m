%%  getRefinementPars %%

%% Return structure of parameters to pass to refineOrientation
%%inputs:
%%  model        - reconstruction volume to use as reference
%%  projections  - N x N x num_projections array of input projections
%%  euler_angles - num_projections x 3 array of phi,theta,psi Euler angles 
%%                  for each projection
%%  varargin     - additional parameters in 'key','value' pairs. A list of options can be found below in %% Options

%%outputs:
%%  refinement_pars - structure containing parameters ready for refineOrientation

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function refinement_pars = getRefinementPars(model,projections,euler_angles,varargin)


if mod(length(varargin),2) ~= 0
    error('Additional argument list not divisible by 2. Options should be ''key'',''value'' pairs.')
end

% default values for options
ang_step           = 1;
ang_range          = 5;
bin_factor         = 1;
compare_func       = @alignByRFactor;
oversampling_ratio = 2;
num_refinements    = 1;
model_filename   = 'data/vesicle.mat';
pj_filename      = 'data/projections.mat';
angle_filename   = 'data/angles.mat';
support_filename = 'data/support.mat';
recon_filename   = 'results/GENFIRE_rec.mat';

%% Options. Set to default values then override with any user-provided parameters
refinement_pars.model               = model; %reconstruction volume to use as reference
refinement_pars.projections         = projections; %projections used to produce model
refinement_pars.euler_angles        = euler_angles;%phi,theta,psi angles
refinement_pars.ang_search_range    = -ang_range:ang_step:ang_range;%vector of angular displacements to search
refinement_pars.phi_search_range    = refinement_pars.ang_search_range;%vector of angular displacements to search phi
refinement_pars.theta_search_range  = refinement_pars.ang_search_range;%vector of angular displacements to search theta
refinement_pars.psi_search_range    = refinement_pars.ang_search_range;%vector of angular displacements to search psi
refinement_pars.bin_factor          = bin_factor;%integer factor by which to bin the model and projections
refinement_pars.compare_func        = compare_func;%comparison function for input projection and 
% calculated backprojection. Must be of the form [metric, new_center_x, new_center_y] = function(input_img,calc_img) where
% metric is the value for R-factor, Xcorr, etc and new_center_x is the
% optimal location of the center found for input_img based upon comparison
% with calc_img
refinement_pars.maximize            = false;%determines whether metric from compare_func should be maximized or minimized
refinement_pars.use_parallel        = true;%use parallel (parfor) where applicable
refinement_pars.oversampling_ratio  = oversampling_ratio;%pad zeros to model for projection calculation to match this OR
refinement_pars.model_filename      = model_filename; %filename of model
refinement_pars.pj_filename         = pj_filename; % filename of projections
refinement_pars.angle_filename      = angle_filename; %filename of euler_angles
refinement_pars.support_filename    = support_filename;%filename of support
refinement_pars.recon_filename      = recon_filename;%filename of reconstruction

% Apply user-provided options
par_number = 1;
while par_number < length(varargin)
    if isfield(refinement_pars,varargin{par_number})
        if strcmp(varargin{par_number},'ang_search_range')
            refinement_pars = setfield(refinement_pars,varargin{par_number},varargin{par_number+1});
            refinement_pars = setfield(refinement_pars,'phi_search_range',varargin{par_number+1});
            refinement_pars = setfield(refinement_pars,'theta_search_range',varargin{par_number+1});
            refinement_pars = setfield(refinement_pars,'psi_search_range',varargin{par_number+1});
        else
            refinement_pars = setfield(refinement_pars,varargin{par_number},varargin{par_number+1});
        end
        par_number = par_number + 2;
    else
        error('Invalid option %s provided.',varargin{par_number})
    end
end

if size(model,1) ~= size(projections,1)
    error('Projections and model should have the same array width')
else
    refinement_pars.model            = bin(refinement_pars.model,refinement_pars.bin_factor,3);
    refinement_pars.projections      = bin(refinement_pars.projections,refinement_pars.bin_factor,2);

end
end

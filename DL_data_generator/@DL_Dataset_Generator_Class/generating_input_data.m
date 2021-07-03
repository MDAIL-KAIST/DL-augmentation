%%%
% J.Lee, KAIST (Korea), 2020.
% Y.Yang, Multi-Dimensional Atomic Imaging Lab, KAIST
% Generating input data 

function obj = generating_input_data(obj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% input parameter %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% # of input file (~10000) 
Number_data_start = obj.Number_data_start;
Number_data_final = obj.Number_data_final;
Number_data=Number_data_final-Number_data_start+1; % total Number

Res = obj.Res;
volsize = obj.volsize;

TARGET_PATH = obj.TARGET_PATH;
TARGET_prefix = obj.TARGET_prefix;
INPUT_PATH = obj.INPUT_PATH;
INPUT_prefix = obj.INPUT_prefix;


% angle parameter for perfect 3D tomogram from atomic structure
input_Random_missing_wedge = obj.input_Random_missing_wedge;
input_theta_angle_range = obj.input_theta_angle_range;

% position to 3d vol (input parameter)
Heights = obj.input_Height;                 % Atom intensity heights
AtomNumbers = obj.input_Atomnumbers;        % Atomic numbers for Fe and Pt
CropHalfWidth = obj.input_CropHalfWidth;    % HalfCropWidth for inserting Gaussian-shaped atom into the volume. 
                                            % This parameter should be large enough so that the atoms will not be sharp-cropped, 
                                            % This para will depend on the pixel size and B factors.
Bfactors = obj.input_Bfactors;              % Gaussian Widths (Debye-Waller factor)            
Bfactors_sigma = obj.input_Bfactors_sigma;  % random Bfactor
Poisson_noise = obj.input_Poisson_noise;    % Adding Gaussian+Poisson noise % yes:1 no:0

Angle_noise_option = obj.input_Angle_noise_option; % Adding Proj angle noise % 0: off, 1: theta 2: all
angle_noise_range = obj.input_Angle_noise_range;         % 0.5 degree uniformly random

%
Vol_normalize_YN = obj.input_Vol_normalize_YN;
Vol_normalize_atomtype = obj.input_Vol_normalize_atomtype;
                    
% 3d vol to 2d images (input parameter)
phi=0;
theta=0;
psi=0;

vec1=[0 0 1];
vec2=[0 1 0];
vec3=[1 0 0];

% GF (input parameter)
% using gpu
rng('shuffle');

%%% calculating fa (electron scattering factor in 3D fourier space)
fa = Func_generate_fa(AtomNumbers,volsize,Res);

%%% preventing repeated tmp files. using time.
dt = datestr(now,'dd_mm_yyyy-HH.MM.SS');

%%

for i=Number_data_start:1:Number_data_final
    %%%% position to 3d vol
    model = importdata(sprintf("%s/%s_coordinate_%d.mat",TARGET_PATH,TARGET_prefix,i));
    atomtype = importdata(sprintf("%s/%s_atomtype_%d.mat",TARGET_PATH,TARGET_prefix,i));
    
    atomlist = find(atomtype~=0);
    model = model(:,atomlist);
    atomtype = atomtype(atomlist);
    
    tic;
    % making vol
    if sum(abs(Bfactors_sigma(:)))~=0
        Vol = single(My_create_vol_from_model_rb_efficient(model, Func_atomnumber2atomtype(atomtype), Heights, Bfactors, Bfactors_sigma, AtomNumbers, volsize, Res, CropHalfWidth ,fa));
    else
        Vol = single(My_create_vol_from_model_exact_HB_efficient(model, Func_atomnumber2atomtype(atomtype), Heights, Bfactors, AtomNumbers, volsize, Res, CropHalfWidth, fa));
    end
    GenTime = toc;
    GenTime = round(10*GenTime)./10;
    fprintf('Generating ideal tomogram has been completed in %.12g seconds.\n\n',GenTime);

    
    %%%% generating randonmly anglelist
    raw_lb_angle = min(input_theta_angle_range);
    raw_ub_angle = max(input_theta_angle_range);
    angle_step = (max(input_theta_angle_range)-min(input_theta_angle_range))/(length(input_theta_angle_range)-1);
    
    if input_Random_missing_wedge~=0
        lb_angle=raw_lb_angle+angle_step*randi([1, round(input_Random_missing_wedge/angle_step)]); % lower bound angle
        ub_angle=raw_ub_angle-angle_step*randi([1, round(input_Random_missing_wedge/angle_step)]); % upper bound angle
        anglelist=lb_angle:angle_step:ub_angle;
    else
        anglelist=input_theta_angle_range;
    end
    
    tic;
    %%%% 3d vol to 2d images
    Car_Projs = zeros(size(Vol,1),size(Vol,2),length(anglelist));
    
    %Vol_k=my_fft(Vol);
    j=1;
    for angle = anglelist
        theta = angle; % input
        %[Car_Projs(:,:,j),~] = calculate3Dprojection_interp_general(Vol_k,phi,theta,psi,vec1,vec2,vec3);
        Car_Projs(:,:,j) = Func_calculate3Dprojection_car(Vol,phi,theta,psi,vec1,vec2,vec3);
        j = j+1;
    end
    Car_Projs(Car_Projs<0) = 0;
    %clear Vol_k
    GenTime = toc;
    GenTime = round(10*GenTime)./10;
    fprintf('Generating ideal tilt series has been completed in %.12g seconds.\n\n',GenTime);

    %%% adding Poisson noise %%%
    if Poisson_noise == 1
        %Car_Projs=poissrnd(Car_Projs)+normrnd(mu,sigma,size(Car_Projs));
        Car_Projs = poissrnd(Car_Projs);
        Car_Projs(Car_Projs<0) = 0;
        fprintf("adding Possion noise \n");
    end
    %%%
    
    %%% adding angle noise %%%
    if Angle_noise_option~=0
        if Angle_noise_option==1 
            Angles = zeros(length(anglelist),3);
            anglelist = anglelist + angle_noise_range*(rand(size(anglelist))-1/2);
            Angles(:,2) = anglelist;
            fprintf("adding angle noise (theta) \n");
        elseif Angle_noise_option==2
            Angles = zeros(length(anglelist),3);
            Angles(:,2) = anglelist;
            Angles = Angles + angle_noise_range*(rand(size(Angles))-1/2); 
            fprintf("adding angle noise (all) \n");
        end
    else
        Angles = zeros(length(anglelist),3);
        Angles(:,2) = anglelist;
    end
    %%%
    
    input_info.Angles = Angles;
    save(sprintf("%s/%s_info_%d.mat",INPUT_PATH,INPUT_prefix,i),'input_info','-v6');
    save(sprintf('./tmp_projs_%s.mat',dt),'Car_Projs');
    save(sprintf('./tmp_angles_%s.mat',dt),'Angles');
    
    %%%%% 2d -> 3d vol (GF)
    % input parameter
    pj_filename = sprintf('./tmp_projs_%s.mat',dt);
    angle_filename = sprintf('./tmp_angles_%s.mat',dt);
    %prefixstr = 'GF_test_dataset_result_Car_IR01';
    
    
    % GENFIRE start (Collecting input parameter)
    GENFIRE = GENFIRE_Reconstructor();

    GENFIRE.filename_Projections = pj_filename;
    GENFIRE.filename_Angles = angle_filename ;
    %GENFIRE.filename_Support = './data/Emptysupport.mat'; 
    %GENFIRE.filename_Results = sprintf('%s.mat',prefixstr);

    GENFIRE.numIterations = obj.GENFIRE.numIterations; % number of iterations
    %GENFIRE.pixelSize = obj.GENFIRE.pixelSize;   % pixel size (in Angstrom)
    GENFIRE.oversamplingRatio = obj.GENFIRE.oversamplingRatio; 
    GENFIRE.griddingMethod = obj.GENFIRE.griddingMethod;  % 1 for FFT, 2 for DFT
    GENFIRE.allowMultipleGridMatches = obj.GENFIRE.allowMultipleGridMatches; % only for FFT, this parameter will be ignored for DFT interpolation
    GENFIRE.constraintEnforcementMode = obj.GENFIRE.constraintEnforcementMode; % Determines the Fourier expansion / supression behavior, "3" turns it off.
    GENFIRE.interpolationCutoffDistance = obj.GENFIRE.interpolationCutoffDistance; % interpolation radius 
    GENFIRE.constraintPositivity = obj.GENFIRE.constraintPositivity; 
    GENFIRE.constraintSupport = obj.GENFIRE.constraintSupport;
    %GENFIRE.ComputeFourierShellCorrelation = obj.GENFIRE.ComputeFourierShellCorrelation; % FSC calculation flag
    %GENFIRE.numBins = obj.GENFIRE.numBins; % number of bins for FSC claculation
    %GENFIRE.percentValuesForRfree = obj.GENFIRE.percentValuesForRfree; % parameter for Rfree calculation
    %GENFIRE.numBinsRfree = obj.GENFIRE.numBinsRfree;
    %GENFIRE.doCTFcorrection = obj.GENFIRE.doCTFcorrection;
    %GENFIRE.CTFThrowOutThreshhold = obj.GENFIRE.CTFThrowOutThreshhold;
    %GENFIRE.calculate_Rfree = obj.GENFIRE.calculate_Rfree; % flag for calculation of Rfree
    GENFIRE.DFT_doGPU = obj.GENFIRE.DFT_doGPU; % flag for using GPU for gridding
    %GENFIRE.REC_doGPU = obj.GENFIRE.REC_doGPU; % flag for using GPU for reconstruction
    %GENFIRE.particleWindowSize = [];
    %GENFIRE.phaseErrorSigmaTolerance = [];
    %GENFIRE.vector1 = [0 0 1];
    %GENFIRE.vector2 = [0 1 0];
    GENFIRE.vector3 = [1 0 0];


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
    
    %run reconstruction
    GENFIRE = reconstruct(GENFIRE);
    
    GF_Vol = GENFIRE.reconstruction;
    
    %%% normalize
    if Vol_normalize_YN==1
        Padded_Vol=zeros(volsize+20,volsize+20,volsize+20);
        Padded_Vol(11:volsize+10,11:volsize+10,11:volsize+10)=GF_Vol;
        
        % extracting only element for normalizing
        normalize_atomtype_index=Func_atomnumber2atomtype(atomtype);
        normalize_atomtype_index=normalize_atomtype_index==Vol_normalize_atomtype;
        
        normalize_input_coordinates = model./Res+round((volsize+1)/2) + 10;
        intensity_3x3 = zeros(1,length(normalize_input_coordinates));

        for j=1:size(normalize_input_coordinates,2)
            curr_pos = round(normalize_input_coordinates(:,j));
            box_3x3 = Padded_Vol((curr_pos(1)-1):(curr_pos(1)+1),(curr_pos(2)-1):(curr_pos(2)+1),(curr_pos(3)-1):(curr_pos(3)+1));
            intensity_3x3(j) = sum(box_3x3(:));
        end

        normalize_factor = nanmean(intensity_3x3)/(3^3);
        clear Padded_Vol

        % normalized volume for input data
        GF_Vol = single(GF_Vol./normalize_factor);
    end
    
    save(sprintf("%s/%s_%d.mat",INPUT_PATH,INPUT_prefix,i),'GF_Vol','-v6');
    clear GENFIRE
    
    fprintf("%d th is finished. \n",i);
       
end

delete(sprintf('./tmp_projs_%s.mat',dt));
delete(sprintf('./tmp_angles_%s.mat',dt));

fprintf(sprintf("saving folder: %s, filename: %s \n",INPUT_PATH,INPUT_prefix));

fprintf("%d (%d-%d) input data have been generated \n",Number_data,Number_data_start,Number_data_final);

end

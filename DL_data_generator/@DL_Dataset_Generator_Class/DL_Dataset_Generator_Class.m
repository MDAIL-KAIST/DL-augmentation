% J.Lee, KAIST, 2021.
% Generate training/test dataset for deep learning augmentation.

classdef DL_Dataset_Generator_Class

    properties
        
        %%% common parameters
        Res;
        volsize;
        
        % number of data
        Number_data_start = 1;
        Number_data_final = 1;
        %Random_Seed;
        Atomnumbers;
        
        %
        INPUT_prefix='input';
        TARGET_prefix='target';
        INPUT_PATH;
        TARGET_PATH;
        
        %%% target data generator values
        target_gen_YN = 0;
        
        % 
        target_cell_vector_a;
        target_cell_vector_b;
        target_cell_vector_c;
        target_init_atoms;
        target_init_basis;
        
        % random effect
        target_Random_shape = 1; % on/off
        target_Random_defect = 0; % percent
        target_Random_rotation = 180; % degree
        target_Random_shift = 22; % pm
        target_Random_traslation = 1; % pixel
        target_Random_shuffling_atomtype = 0;
        
        % parameter for amorphous structure
        target_atomdensity = 4/3.912^3 ; % target atomdensity
        target_minimum_distance = 2.0 ; % minimum distance between atoms (Angstron)
        target_getting_percent = 0.99; % generating atoms until getting target percent of atomdensity

        % target Gaus vol
        target_gen_Gaus_YN = 1; 
        target_Height;
        target_Bfactors;
        target_Atomnumbers;
        target_CropHalfWidth = 4;
        
 
        %%% input data generator values
        input_data_YN = 0;
        
        % generating 3D tomogram from atomic model
        input_Random_missing_wedge
        input_theta_angle_range
        input_Angle_noise_option
        input_Angle_noise_range
        
        input_Height;
        input_Bfactors;
        input_Bfactors_sigma;
        input_Atomnumbers;
        input_CropHalfWidth = 8;        
        input_Poisson_noise = 1;
        input_Vol_normalize_YN = 0;
        input_Vol_normalize_atomtype = 1;
        
        % projection_option
        input_Projs_method = 1; % 1:fft, 2:real
        input_Projs_GPU = 0; % 1: on, 0: off
        
        % reconstruction(GENFIRE) parameters    
        GENFIRE;
        
    end
 
    methods
        function obj = DL_Dataset_Generator_Class()
            obj.GENFIRE = GENFIRE_Reconstructor();
        end
        
        function MAKEDIR(obj)
            if ~isdir(obj.INPUT_PATH)
                mkdir(obj.INPUT_PATH);
                fprintf('mkdir: %s \n',obj.INPUT_PATH);
            end
            if ~isdir(obj.TARGET_PATH)
                mkdir(obj.TARGET_PATH);
                fprintf('mkdir: %s \n',obj.TARGET_PATH);
            end
        end
        
        function obj = GEN_target_data(obj)
            if obj.target_gen_YN==1
                obj=generating_target_data(obj);
            end
        end
        
        function obj = GEN_target_data_amorphous(obj)
            if obj.target_gen_YN==1
                obj=generating_target_data_amorphous(obj);
            end
        end
        
        function obj = GEN_input_data(obj)
            if obj.input_data_YN==1
                obj=generating_input_data(obj);
            end
        end
        
        % test
        
        function obj = GEN_target_data_new(obj)
            if obj.target_gen_YN==1
                obj=generating_target_data_v1(obj);
            end
        end
        
        function obj = GEN_input_data_new(obj)
            if obj.input_data_YN==1
                obj=generating_input_data_v1(obj);
            end
        end        
        

        
    end
    
end
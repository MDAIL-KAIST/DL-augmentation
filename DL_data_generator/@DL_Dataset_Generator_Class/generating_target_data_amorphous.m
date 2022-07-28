%%%
% J.Lee, KAIST (Korea), 2020.
% Y.Yang, Multi-Dimensional Atomic Imaging Lab, KAIST
% Generating target data (for an amorphous structure)

function obj=generating_target_data_amorphous(obj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% input parameter %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% # of input file (~10000) 
Number_data_start = obj.Number_data_start;
Number_data_final = obj.Number_data_final;
Number_data=Number_data_final-Number_data_start+1; % total Number

TARGET_PATH = obj.TARGET_PATH;
TARGET_prefix = obj.TARGET_prefix;

Res = obj.Res;
volsize = obj.volsize;

%%% option for translation, noise, shifting each atom, rotation ...%%%%
Random_shape = obj.target_Random_shape; % extracting random shape , yes:1  no:0; 
atomdensity = obj.target_atomdensity; % target density
minimum_distance = obj.target_minimum_distance; % minimum distance between atoms
getting_percent = obj.target_getting_percent; % generating atoms until getting target percent of atomdensity

%%% Generating target data
gen_Gaus_YN = obj.target_gen_Gaus_YN; % yes:1  no:0;

% gaussian (input parameter)
Heights = obj.target_Height;              % Atom intensity heights (e.g. Fe and Pt -> [26 78]); 
Bfactors = obj.target_Bfactors;           % Gaussian Widths (Debye-Waller factor); 
AtomNumbers = obj.target_Atomnumbers;     % Atomic numbers (e.g. Fe and Pt -> [26 78]);
CropHalfWidth = obj.target_CropHalfWidth; % HalfCropWidth for inserting Gaussian-shaped atom into the volume.  %4
                                          % This para will depend on the pixel size and B factors.

rng('shuffle');
% random shape (input parameter)
volumeMinMax=[volsize.^3/6 volsize.^3/4]; % 8, 6
sphereRadiusMinMax=[round(volsize/8) round(volsize/4)]; %12 4
smoothingKernelSize=9; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for i=Number_data_start:1:Number_data_final
    %%% generating amorphous model
    atomcoordinates = Func_generating_pos_3D_rand(volsize, Res, minimum_distance, atomdensity, getting_percent, 1);
    atomtype = zeros(1,length(atomcoordinates))+1;
    %
    atomcoordinates_res = atomcoordinates./Res;
    OutInd = find(atomcoordinates_res(1,:) < 1 | atomcoordinates_res(1,:) > volsize | ...
        atomcoordinates_res(2,:) < 1 | atomcoordinates_res(2,:) > volsize | ...
        atomcoordinates_res(3,:) < 1 | atomcoordinates_res(3,:) > volsize);
    atomcoordinates_res(:,OutInd) = [];
    atomtype(:,OutInd) = [];
    atomcoordinates_round=round(atomcoordinates_res);
    atomcoordinates=Res*(atomcoordinates_res-round((volsize+1)/2));
    
    
    % generating hole randomly
    new_rand_coordinate=atomcoordinates;
    new_rand_atomtype=atomtype;
    
    % extracting random shape
    if Random_shape==1
        RandVol = My_generate_random_volume(volsize,volumeMinMax, sphereRadiusMinMax, smoothingKernelSize);
  
        Ind = sub2ind(size(RandVol),atomcoordinates_round(1,:),atomcoordinates_round(2,:),atomcoordinates_round(3,:));
        GoodAtomInd = find(RandVol(Ind)>0);
        
        new_rand_coordinate = new_rand_coordinate(:,GoodAtomInd);
        new_rand_atomtype = new_rand_atomtype(GoodAtomInd);
        
    else
        Norms = sqrt(sum(f_pos.^2,1));
        ParticleRadius=volsize*Res/2*0.8; % default value
        ParticleInd  = find(Norms < ParticleRadius);
        new_rand_coordinate = f_pos(:,ParticleInd);
        new_rand_atomtype = f_att(ParticleInd);
    end     
    

    atomlist=find(new_rand_atomtype~=0);
    new_rand_coordinate=new_rand_coordinate(:,atomlist);
    new_rand_atomtype=new_rand_atomtype(atomlist);
    
        
    if gen_Gaus_YN==1
        target = Func_create_voxel_gau(new_rand_coordinate, Func_atomnumber2atomtype(new_rand_atomtype), Heights, Bfactors, AtomNumbers, volsize, Res,CropHalfWidth);
        target = single(target);
    end
    
    %    
    save(sprintf("%s/%s_%d",TARGET_PATH,TARGET_prefix,i),"target",'-v6');
    save(sprintf("%s/%s_coordinate_%d",TARGET_PATH,TARGET_prefix,i),"new_rand_coordinate",'-v6');
    save(sprintf("%s/%s_atomtype_%d",TARGET_PATH,TARGET_prefix,i),"new_rand_atomtype",'-v6');
    
    fprintf("%d th is finished. \n",i);
end

fprintf(sprintf("saving folder: %s, filename: %s \n",TARGET_PATH,TARGET_prefix));

fprintf("%d (%d-%d) target data have been generated \n",Number_data,Number_data_start,Number_data_final);

end


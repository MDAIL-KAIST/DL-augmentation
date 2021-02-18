%%%
% J.Lee, KAIST (Korea), 2020.
% Y.Yang, Multi-Dimensional Atomic Imaging Lab, KAIST
% Generating target data 

function obj=generating_target_data(obj)

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
Random_defect = obj.target_Random_defect; % defect, percent
Random_angle = obj.target_Random_rotation; % rotation 
Random_traslation = obj.target_Random_traslation; % random translation, yes:1  no:0;
Random_shift = obj.target_Random_shift; % shifting each atom, pm scale
Random_shuffling_atomtype = obj.target_Random_shuffling_atomtype; % shuffling atom randomly

%%% cell_vector + basis
cell_vector_a = obj.target_cell_vector_a;
cell_vector_b = obj.target_cell_vector_b;
cell_vector_c = obj.target_cell_vector_c;

init_att = obj.target_init_atoms; % atomtype
init_pos = obj.target_init_basis; % basis


%%% Generating target data
gen_Gaus_YN = obj.target_gen_Gaus_YN; % yes:1  no:0;

% gaussian (input parameter)
Heights = obj.target_Height;              % Atom intensity heights (e.g. Fe and Pt -> [26 78]); 
Bfactors = obj.target_Bfactors;           % Gaussian Widths (Debye-Waller factor); 
AtomNumbers = obj.target_Atomnumbers;     % Atomic numbers (e.g. Fe and Pt -> [26 78]);
CropHalfWidth = obj.target_CropHalfWidth; % HalfCropWidth for inserting Gaussian-shaped atom into the volume.  %4
                                          % This para will depend on the pixel size and B factors.

rng('shuffle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

lattice_para=min([norm(cell_vector_a) norm(cell_vector_b) norm(cell_vector_c)]);
numatom_per_UC = length(init_att);

% cell_vector matrix
cell_vector=[cell_vector_a; cell_vector_b; cell_vector_c]';

%%

% number of unit cells (UCs) to initially simulate
Natoms_per_dim = round(volsize * Res * 0.6 / lattice_para);

% initialize position and atomtype
f_pos = repmat(init_pos,[1 (Natoms_per_dim*2+1)^3]);
f_att = repmat(init_att,[1 (Natoms_per_dim*2+1)^3]);

X = -1*Natoms_per_dim:1:Natoms_per_dim;
[X,Y,Z] = meshgrid(X,X,X);

% set proper lattice vectors for all UC points
for i=1:numatom_per_UC
    curr_ind = (1:numatom_per_UC:(Natoms_per_dim*2+1)^3*numatom_per_UC) + (i-1);
    f_pos(:,curr_ind) = f_pos(:,curr_ind) + [X(:)'; Y(:)'; Z(:)'];
end

% apply lattice parameter to create atomic coordinates from the UC
% structure
f_pos = cell_vector*f_pos;

atomcoordinates = f_pos;
atomtype = f_att;


%% 
atomcoordinates_res = atomcoordinates/Res + round((volsize+1)/2);
OutInd = find(atomcoordinates_res(1,:) < 1 | atomcoordinates_res(1,:) > volsize | ...
    atomcoordinates_res(2,:) < 1 | atomcoordinates_res(2,:) > volsize | ...
    atomcoordinates_res(3,:) < 1 | atomcoordinates_res(3,:) > volsize);
atomcoordinates_res(:,OutInd) = [];
atomtype(:,OutInd) = [];
atomcoordinates_round=round(atomcoordinates_res);

atomcoordinates=Res*(atomcoordinates_res-round((volsize+1)/2));

                    
% random shape (input parameter)
volumeMinMax=[volsize.^3/8 volsize.^3/5];
sphereRadiusMinMax=[round(volsize/8) round(volsize/4)]; %12 4
smoothingKernelSize=9; 

%%
for i=Number_data_start:1:Number_data_final
    % generating hole randomly
    new_rand_coordinate=atomcoordinates;
    new_rand_atomtype=atomtype;
    
    % extracting random shape
    if Random_shape==1
        RandVol = My_generate_random_volume(volsize,volumeMinMax, sphereRadiusMinMax, smoothingKernelSize);
  
        Ind = sub2ind(size(RandVol),atomcoordinates_round(1,:),atomcoordinates_round(2,:),atomcoordinates_round(3,:));
        GoodAtomInd = find(RandVol(Ind)>0);
        
        new_rand_coordinate = new_rand_coordinate(:,GoodAtomInd);
        new_rand_atomtype = new_rand_atomtype(:,GoodAtomInd);
        
    else
        Norms = sqrt(sum(f_pos.^2,1));
        ParticleRadius=volsize*Res/2*0.8; % default value
        ParticleInd  = find(Norms < ParticleRadius);
        new_rand_coordinate = f_pos(:,ParticleInd);
        new_rand_atomtype = f_att(ParticleInd);
    end     
    
    % defining # of atom
    Nofatom=length(new_rand_atomtype);
    
    % shift each atom (rms_set) pixel (sub-pixel)
    if Random_angle~=0
        direction_rand_theta=pi/2*rand(1);
        direction_rand_phi=2*pi*rand(1);
        direction_vec=[cos(direction_rand_theta)*cos(direction_rand_phi),...
            cos(direction_rand_theta)*sin(direction_rand_phi), sin(direction_rand_theta)]; % random_direction vector
        
        random_theta=Random_angle*(2*rand(1)-1); % putting random angle
        Rot_matrix=MatrixQuaternionRot(direction_vec,random_theta);
        
        new_rand_coordinate=Rot_matrix*new_rand_coordinate;
        
        %
        target_info.Random_rotation_direction=direction_vec;
        target_info.Random_rotation_angle=random_theta;
    else
        target_info.Random_rotation_direction=[];
        target_info.Random_rotation_angle=[];
    end
    
    % shift each atom (rms_set) pixel (sub-pixel)
    if Random_shift~=0
        endflag=0;
        while ~endflag
           shift_para=2*Random_shift*rand(3,Nofatom)-Random_shift; 
           mean_shift=mean((sum((shift_para).^2)).^0.5);
           if mean_shift>Random_shift*0.9 && mean_shift<Random_shift*1.1
               endflag=1;
           end    
        end
        new_rand_coordinate=new_rand_coordinate+shift_para;
    end
    
    % translation every atom xx pixel
    if Random_traslation~=0
        random_tr=Random_traslation*rand(3,1)*Res;
        new_rand_coordinate=new_rand_coordinate+random_tr;
    end
    
    
    % making hole = removing atom randomly
    if Random_defect~=0
        if round(Nofatom*Random_defect)>=1
            Nofhole=randi([1 round(Nofatom*Random_defect/100)],1,1); % making hole (# of holes) 
            for j=1:1:Nofhole
                changing_c=randi([1 Nofatom+1-j],1,1);
                new_rand_coordinate(:,changing_c)=[0;0;0];
                new_rand_atomtype(changing_c)=0; 
            end
        end
    end
 

    atomlist=find(new_rand_atomtype~=0);
    new_rand_coordinate=new_rand_coordinate(:,atomlist);
    new_rand_atomtype=atomtype(atomlist);
    
    if Random_shuffling_atomtype==1
        new_rand_atomtype=new_rand_atomtype(randperm(length(new_rand_atomtype)));
    end
        
    if gen_Gaus_YN==1
        target = Func_create_voxel_gau(new_rand_coordinate, Func_atomnumber2atomtype(new_rand_atomtype), Heights, Bfactors, AtomNumbers, volsize, Res,CropHalfWidth);
        target = single(target);
    end
    
    %    
    save(sprintf("%s/%s_%d",TARGET_PATH,TARGET_prefix,i),"target",'-v6');
    save(sprintf("%s/%s_coordinate_%d",TARGET_PATH,TARGET_prefix,i),"new_rand_coordinate",'-v6');
    save(sprintf("%s/%s_atomtype_%d",TARGET_PATH,TARGET_prefix,i),"new_rand_atomtype",'-v6');
    save(sprintf("%s/%s_info_%d",TARGET_PATH,TARGET_prefix,i),"target_info",'-v6');
    
    fprintf("%d th is finished. \n",i);
end

fprintf(sprintf("saving folder: %s, filename: %s \n",TARGET_PATH,TARGET_prefix));

fprintf("%d (%d-%d) target data have been generated \n",Number_data,Number_data_start,Number_data_final);

end


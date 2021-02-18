clear all
clc
%%
% J.Lee, KAIST, 2021.
% Y.Yang, Multi-Dimensional Atomic Imaging Lab
% Pre-proccesing for DL augmentation input data
% Apply zero-padding (match DL augmentation volume size)
% and real interpolation (upsampling(3x3x3 voxel) + binning(average)

%%
%Vol_INPUT_PATH = '~/Dropbox/share/mdail_test/data/exp_data/Pt_inputdata/calib_result';
Vol_INPUT_PATH = './raw_reconstructure/';
Vol_OUTPUT_PATH = './Preproccesing_tomogram_for_DL_aug/';
ESTvol=importdata(sprintf('%s/Pt_input_1.mat',Vol_INPUT_PATH));

%%
volsize=144;
ESTvol=single(My_paddzero(ESTvol,[volsize,volsize,volsize]));

%
OV=3;
ESTvol = Func_real_intepolation_OV_sampling(ESTvol,OV);

%
save(sprintf('%s/Pt_input_1_real_intepolation_zero_padding.mat',Vol_OUTPUT_PATH),'ESTvol');

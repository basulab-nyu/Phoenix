clear;
%% load all modules
addpath(genpath(pwd))
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/Phoenix')) %path to CNMF 

%% load files 
% Folder with stack images of diferent sessions
foldername='/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/M1/FOV1/'; %path to image folder
cd(foldername);
listfiles = dir('*.tif'); %will look for all the tif files
% !!! Files are listed in alphabetical order !!!
for i= 1:length(listfiles)
files{i}=fullfile(foldername, listfiles(i).name);
end
% Import ROI
FOV = [512,512]; % Image resolution 
d1=FOV(1);
d2=FOV(2);
ROI_file='/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/M1/ROI/M1_FOV1.zip'; %path to ROI file (.zip)
% Need ReadImageJROI.m script
[a,ROI] = ReadImageJROI(ROI_file,[d1,d2]);
% Create Structure
input.foldername=foldername;
input.image=files;
input.list=listfiles;
input.ROI.file=ROI;
input.ROI.a=a;
input.param.FOV=FOV;

%% Set parameters soma
K = 1;                                           % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron)
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
options = CNMFSetParms(...
                       'd1',d1,'d2',d2,...                        % dimensions of datasets
                       'search_method','dilate','dist',3,...       % search locations when updating spatial components
                       'deconv_method','constrained_foopsi',...    % activity deconvolution method
                       'temporal_iter',2,...                       % number of block-coordinate descent steps
                       'fudge_factor',0.98,...                     % bias correction for AR coefficients
                       'merge_thr',merge_thr,...                    % merging threshold
                       'gSig',tau...
                       );
% Create Structure                
input.param.options=options;
input.param.tau=tau;
input.param.p=p;
input.param.merge_thr=merge_thr;
                  
%% CNMF
input.refine=0; %refine components 
%Session to process
%session=1; 
%OR
%all session
session=1:size(files,2);

for i=1:length(session)
[output]= CNMF_noparpool(input, session(i)); 
end






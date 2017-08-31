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

%% Set parameters dendrites
K = 2;                                           % number of components to be found
tau = [];                                         % std of gaussian kernel (size of neuron - not needed for dendritic data) 
p = 0;                                            % No AR dynamics for dendritic data
merge_thr = 0.8;                                  % merging threshold
options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                           % dimensions of datasets
    'init_method','HALS',...                      % initialize algorithm with plain NMF  
    'max_iter_hals_in',50,...                     % maximum number of iterations
    'search_method','dilate',...                  % search locations when updating spatial components
    'temporal_iter',2,...                         % number of block-coordinate descent steps 
    'merge_thr',0.9,...                           % merging threshold
    'conn_comp',false,...                         % do not limit to largest connected component for each found component
    'maxthr',0.05...                              % for every component set pixels that are below max_thr*max_value to 0 
    );
% Create Structure                
input.param.options=options;
input.param.tau=tau;
input.param.p=p;
input.param.merge_thr=merge_thr;
                  
%% CNMF
input.refine=0; %refine components 
input.refine=0; %refine components 
%Session to process
%session=1; 
%OR
%all session
session=1:size(files,2);

for i=1:length(session)
[output]= CNMF_noparpool(input, session(i)); 
end






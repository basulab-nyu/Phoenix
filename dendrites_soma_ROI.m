clear;
gcp;
addpath(genpath(pwd))
addpath(genpath('/Users/martial/GitHub/Phoenix'))  

%% load file 
nam = '/Users/martial/Documents/test/1_M991_CTRL_MC.tif'; %Images to open
%Import ROI
FOV = [512,512]; % Image resolution 
d1=FOV(1);
d2=FOV(2);
dend='/Users/martial/Desktop/PSAM_ROI/M991/soma&dendrites/1_2_dendrites.zip'; %path to ROI file (.zip)
%soma='/Users/martial/Desktop/PSAM_ROI/M991/soma_10micron_polygon.zip'; %path to ROI file (.zip)
[~,ROI_dend] = ReadImageJROI(dend,[d1,d2]);
[~,ROI_soma] = ReadImageJROI(soma,[d1,d2]);
%Open image
sframe=1;	% first frame to read (start of the session)
num2read=1000;	% how many frames to read   (lenght of the session)
Y = bigread2(nam,sframe,num2read);
%Y = bigread2(nam);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;  
%% Save structures
input.name=nam;
input.dend=dend;
input.soma=soma;
input.ROI_dend=ROI_dend;
input.ROI_soma=ROI_soma;
%% Data pre-processing
[P,Y] = preprocess_data(Y,0); %set at 0 for dendrites
Cn_max =  correlation_image_max(Y);
%% Start CNMF for dendrites and soma
input.update_dendrite=0 %Update spatial components for dendrites?
CNMF_dendrites(Y,P,Cn_max,input); 
CNMF_soma(Y,P,Cn_max,input);  

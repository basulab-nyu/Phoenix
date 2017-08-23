clear;
gcp;
clear;
addpath(genpath(pwd))
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/ca_source_extraction')) %CNMF 
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/NoRMCorre')) %NorMCorre
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/scripts')) 


%Folder with stack images of diferent sessions 
foldername='/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/M2/FOV1';
cd(foldername);
listhdf5 = dir('*.hdf5'); %will look for all the hdf5 files
% !!! Files are listed in alphabetical order !!!
for i= 1:length(listhdf5)
files{i}=fullfile(foldername, listhdf5(i).name);
end
FOV = [512,512];

%% Set parameters
options_rigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',200,'max_shift',30,'us_fac',50, 'phase_flag', true);
tsub=100; % temporal downsampling factor 

%% MC and Downsample for drawing ROI 
%Create a blank template for the 1st session to MC
template=[];
%Or load a template from previous sessions:
%temp_sess=1; %template from session
%load(([files{temp_sess}(1:end-5),'_MC_output']), 'template');

for i=1:size(files,2) %MC all files
%OR MC only sessions X, Y, ..
%
%sessions=10:12;
%for i=sessions


name=files{i};

% Open stack
disp(['Opening : ' listhdf5(i).name])
Y=read_file(name);
Y=single(Y);
disp(['Done opening : ' listhdf5(i).name])

% Start MC
disp(['Motion Correction : ' listhdf5(i).name])
tic; [M,shifts,template] = normcorre_batch(Y,options_rigid,template);
toc;
disp(['Done Motion Correction : ' listhdf5(i).name])

% Save MC images + shifts and options
disp(['Saving : ' listhdf5(i).name])
tic;
save(fullfile([name(1:end-5),'_MC_output']),'shifts', 'options_rigid', 'template');
file_name=[name(1:end-5),'_MC.tif'];
saveastiff(int16(M), file_name);
disp(['Done saving : ' listhdf5(i).name])
toc;

% Perform temporal downsampling (mean)
ndimsY = ndims(Y)-1;
sY = size(Y);
ds = sY(1:ndimsY);
T = sY(end);
Ts = floor(T/tsub);   
Y_ds = squeeze(mean(reshape(Y(:, :, 1:(Ts*tsub)),ds(1), ds(2), tsub, Ts), 3));

%Perform temporal downsampling (correlation);
for ii=1:Ts;
sfr=(tsub*ii)-tsub+1;
endfr=tsub*ii;
fr=[sfr; endfr];
Ytemp=Y(:,:,fr(1):fr(2));
Cntemp=correlation_image(Ytemp);
Cn_Z(:,:,ii)=Cntemp;
end

%Save as tiff
%name of the temporal downsampled (mean) image
TSub_mean_name=[name(1:end-4),'_Tsub_mean.tif']; 
%name of the temporal downsampled (correlation) image
TSub_CN_name=[name(1:end-4),'_Tsub_CN.tif'];
saveastiff(int16(Y_ds), TSub_mean_name);
saveastiff(Cn_Z, TSub_CN_name);
end



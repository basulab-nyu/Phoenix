%% Script to make a downsample images using mean and correlation.
clear;
%load all modules
addpath(genpath(pwd))
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/ca_source_extraction')) %CNMF 
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/NoRMCorre')) %NorMCorre
addpath(genpath('/ifs/data/basulab/MAD/MATLAB/scripts')) 

%% temporal downsampling factor
tsub=100; 

%% load file
% Folder with stack images of diferent sessions
foldername='/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/M1/FOV1';
cd(foldername);
listfiles = dir('*.tif'); %will look for all the tif files
% !!! Files are listed in alphabetical order !!!
for i= 1:length(listfiles)
files{i}=fullfile(foldername, listfiles(i).name);
end

%% load image
for i=1:size(files,2)
nam=files{i};
Y = bigread2(nam);

%% Temporal downsampling (mean)
ndimsY = ndims(Y)-1;
sY = size(Y);
ds = sY(1:ndimsY);
T = sY(end);
Ts = floor(T/tsub);   
Y_ds = squeeze(mean(reshape(Y(:, :, 1:(Ts*tsub)),ds(1), ds(2), tsub, Ts), 3));

%% Temporal downsampling (correlation);
for ii=1:Ts;
sfr=(tsub*ii)-tsub+1;
endfr=tsub*ii;
fr=[sfr; endfr];
Ytemp=Y(:,:,fr(1):fr(2));
Cntemp=correlation_image(Ytemp);
Cn_Z(:,:,ii)=Cntemp;
end

%% Save as tiff
TSub_mean_name=[nam(1:end-4),'_Tsub_mean.tif'];
TSub_CN_name=[nam(1:end-4),'_Tsub_CN.tif'];
saveastiff(int16(Y_ds), TSub_mean_name);
saveastiff(Cn_Z, TSub_CN_name);
end








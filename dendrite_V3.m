clear;
addpath(genpath('../CNMF'));
%% load file 
nam = '/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/06_19_2017/M1/M1_STD1_STD2_All_mcRigid.tif'; %Images to open
dirnam='/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/06_19_2017/M1/'; %Dir to save files 
%dend='E:\MAD_DATA\CA2_Amigo_FlexGC6f\04_06_2017\M1\all_dend.zip'; 
soma='/ifs/data/basulab/MAD/DATA/CA3_ThyGC6f/06_19_2017/M1/M1_soma.zip'; 

sframe=1;	% first frame to read (start of the session)
num2read=19999;	% how many frames to read   (lenght of the session)
Y = bigread2(nam,sframe,num2read);


%Y = bigread2(nam);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;  
%% Set parameters dendrites
session='dendrites';
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
    'merge_thr',0.8,...                           % merging threshold
    'conn_comp',false,...                         % do not limit to largest connected component for each found component
    'maxthr',0.05...                              % for every component set pixels that are below max_thr*max_value to 0 
    );
%% Data pre-processing
[P,Y] = preprocess_data(Y,p);
%% fast initialization of spatial components using greedyROI and HALS
Cn_max =  correlation_image_max(Y);
Cn = Cn_max;
%% Import ROI from ImageJ
% Need ReadImageJROI.m script
[a,ROI] = ReadImageJROI(dend,[d1,d2]);
%imshow((Cn_max)/max(max((Cn_max))).*nansum(ROI,3));
Ain= sparse(reshape(ROI,512*512,[]));
for k=1:size(ROI,3)
    ROI_bin{k}=ROI(:,:,k);
    Cin(k,:) = ROI_bin{k}(:).'*reshape(Y,[],T)/nnz(ROI_bin{k});
end
Yr = reshape(Y,d,T);
res = reshape(Y,[],T) - Ain*Cin;
fin = mean(res);
for nmfiter = 1:20
bin = max((res*fin')/(fin*fin'),0);    
fin = max((bin'*bin)\(bin'*res),0);
end
%% update spatial components
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed and dendrites
[C,f,P,S,YrA] = update_temporal_components(Yr,Ain,b,Cin,fin,P,options);
%% refine estimates excluding rejected components
[A2,b2,C2] = update_spatial_components(Yr,C,f,[A,b],P,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,P,options);
%% Extract DF/F
[C_df,~] = extract_DF_F(Yr,A2,C2,P2,options);
% center for ROI 
center= com(A2,d1,d2);
% New df/f extract
alpha=0.05;
[expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA2, A2,C2, b2,f2,alpha); % extract DF/F values (optional)
%% display components
[Coor,json_file] = plot_contours(A2,Cn_max,options,1);
plot_dend_components_GUI(Yr,A2,C2,b2,f2,options);

%% Save
namcdf=strcat(dirnam, session, 'Cdf');
namcenterROI=strcat(dirnam,session, 'ROI');
tic; save(namcdf,'C_df', 'expDffMedZeroed' ); toc;
tic; save(namcenterROI,'Coor' ,'json_file', 'center' ); toc;

clearvars -except Y p P d d1 d2 T Cn Cn_max nam dirnam dend soma

%% Set parameters soma
session='CTRL';
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
%% Import ROI from ImageJ
% Need ReadImageJROI.m script
[a,ROI] = ReadImageJROI(soma,[d1,d2]);
%imshow((Cn_max)/max(max((Cn_max))).*nansum(ROI,3));
Ain= sparse(reshape(ROI,512*512,[]));
for k=1:size(ROI,3)
    ROI_bin{k}=ROI(:,:,k);
    Cin(k,:) = ROI_bin{k}(:).'*reshape(Y,[],T)/nnz(ROI_bin{k});
end
Yr = reshape(Y,d,T);
res = reshape(Y,[],T) - Ain*Cin;
fin = mean(res);
for nmfiter = 1:20
bin = max((res*fin')/(fin*fin'),0);    
fin = max((bin'*bin)\(bin'*res),0);
end
%% update spatial components
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
%% update temporal components
P.p = 2;    
[C,f,P,S,YrA] = update_temporal_components(Yr,Ain,b,Cin,fin,P,options);
%% classify components
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);
%% run GUI for modifying component selection (optional, close twice to save values)
%run_GUI = true;
%if run_GUI
%    Coor = plot_contours(A,Cn,options,1); close;
 %   GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
 %   options = GUIout{2};
 %   keep = GUIout{3};    
%end
%% refine estimates excluding rejected components
[A2,b2,C2] = update_spatial_components(Yr,C,f,[A,b],P,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,P,options);
%% Extract DF/F
[C_df,~] = extract_DF_F(Yr,A2,C2,P2,options);
% center for ROI 
center= com(A2,d1,d2);
% New df/f extract
alpha=0.05;
[expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA2, A2,C2, b2,f2,alpha); 
%% display components
[Coor,json_file] = plot_contours(A2,Cn_max,options,1);
plot_components_GUI(Yr,A2,C2,b2,f2,Cn,options)

%% Save
namcdf=strcat(dirnam, session, 'Cdf');
namcenterROI=strcat(dirnam,session, 'ROI');
tic; save(namcdf,'C_df', 'expDffMedZeroed' ); toc;
tic; save(namcenterROI,'Coor' ,'json_file', 'center', 'ROIvars', 'keep' ); toc;


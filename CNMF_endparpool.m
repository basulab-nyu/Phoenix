function [output]= CNMF_endparpool(input, session);  
%% load structure
foldername = input.foldername;
files = input.image;
listhdf5 = input.list;
ROI = input.ROI.file;
a = input.ROI.a;
FOV = input.param.FOV;
options = input.param.options;
tau = input.param.tau;
p = input.param.p;
merge_thr =input.param.merge_thr;
%Session to process
nam=files{session}; 
%sframe=1;
%num2read=5000;
%Y = bigread2(nam,sframe,num2read);
Y = bigread2(nam);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;

%% Data pre-processing
[P,Y] = preprocess_data(Y,p);
%% fast initialization of spatial components using greedyROI and HALS
Cn_max =  correlation_image_max(Y);
Cn = Cn_max;
%Use ROI imported ROI for Ain
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
poolobj = parpool;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 2;    
[C,f,P,S,YrA] = update_temporal_components(Yr,Ain,b,Cin,fin,P,options);
delete(poolobj);

%% classify components
poolobj = parpool;
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);
% run GUI for modifying component selection (optional, close twice to save values)
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
%New df/f extract
alpha=0.05;
[expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA2, A2,C2, b2,f2,alpha);

%% detrend fluorescence and extract DF/F values
options.df_window = 1000; 
[F_dff,F0] = detrend_df_f(A2,b,C2,f2,YrA2,options);

%% deconvolve data
nNeurons = size(F_dff,1);
C_dec = zeros(size(F_dff));
S = zeros(size(F_dff));
kernels = cell(nNeurons,1);
min_sp = 3;    % find spikes resulting in transients above min_sp x noise level

for i = 1:nNeurons
    [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
end


%% display components
[Coor,json_file] = plot_contours(A2,Cn_max,options,1);
plot_components_GUI(Yr,A2,C2,b2,f2,Cn,options)
%center for ROI
center= com(A2,d1,d2);

%Save
disp(['Saving : ' listhdf5(session).name])
tic;
save(fullfile([nam(1:end-6),'Cdf']),'C_dec','C_df','F_dff','expDffMedZeroed','S');
save(fullfile([nam(1:end-6),'ROI']),'Coor' ,'json_file', 'center', 'ROIvars', 'keep','A2','Cn_max', 'options');
disp(['Done saving : ' listhdf5(session).name])

output.expDffMedZeroed= expDffMedZeroed;
output.C_df= C_df;
output.center=center;
output.Coor=Coor;
output.keep=keep;
output.Cn_max=Cn_max;
output.A2=A2;

delete(poolobj);





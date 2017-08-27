
function CNMF_dendrites(Y,P,Cn_max,input);  
%% Import files
[d1,d2,T] = size(Y); 
update=input.update_dendrite;
nam=input.name;
ROI=input.ROI_dend;
d = d1*d2;  
Cn = Cn_max;
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
    'merge_thr',0.9,...                           % merging threshold
    'conn_comp',false,...                         % do not limit to largest connected component for each found component
    'maxthr',0.05...                              % for every component set pixels that are below max_thr*max_value to 0 
    );
%% Import ROI from ImageJ
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
%% Plot imageJ ROI
%plot_contours(Ain,Cn_max,options,1);
imshow((Cn_max)/max(max((Cn_max))).*nansum(ROI,3));
%% update spatial components
if update==1
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
elseif update==0
A=Ain;
b=bin;
end
%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed and dendrites
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);
%% Plot uptdaded ROI
%plot_contours(A,Cn_max,options,1);
%% refine estimates excluding rejected components
%[A2,b2,C2] = update_spatial_components(Yr,C,f,[A,b],P,options);
%[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,P,options);
%% Plot refined ROI
%plot_contours(A2,Cn_max,options,1);
%% Do not refine estimates for dendrites:
A2=A;
b2=b;
C2=C;
P2=P;
f2=f;
S2=S;
YrA2=YrA;
%% Extract DF/F
[C_df,~] = extract_DF_F(Yr,A,C,P,options);
% center for ROI 
center= com(A,d1,d2);
% New df/f extract
alpha=0.05;
[expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA, A,C, b,f,alpha); % extract DF/F values (optional)
%% detrend fluorescence and extract DF/F values
%options.df_window = 1000; 
%[F_dff,F0] = detrend_df_f(A2,b,C2,f2,YrA2,options);
%% deconvolve data
%nNeurons = size(F_dff,1);
%C_dec = zeros(size(F_dff));
%S = zeros(size(F_dff));
%kernels = cell(nNeurons,1);
%min_sp = 3;    % find spikes resulting in transients above min_sp x noise level
%for i = 1:nNeurons
%    [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
%end
%% display components
%[Coor,json_file] = plot_contours(A2,Cn_max,options,1);
%plot_dend_components_GUI(Yr,A2,C2,b2,f2,options);
%% Save
tic;
save(fullfile([nam(1:end-6),session, '_Cdf']),'C_df','expDffMedZeroed');
save(fullfile([nam(1:end-6),session,'_ROI']),'center','A2','C2','b2','f2','Cn_max', 'options');
toc;

end
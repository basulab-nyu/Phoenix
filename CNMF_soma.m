
function CNMF_soma(Y,P,Cn_max,input);  
%% Import files
[d1,d2,T] = size(Y);  
nam=input.name;
ROI=input.ROI_soma;
d = d1*d2; 
Cn = Cn_max;
%% Set parameters soma
session='soma';
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
imagesc((Cn_max)/max(max((Cn_max))).*nansum(ROI,3));

%% update spatial components
%if update==1
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
%elseif update==0
%A=Ain;
%b=bin;
%end

%% update temporal components
P.p = 2;    % set AR to 2
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);
%% Plot uptdaded ROI
%plot_contours(A,Cn_max,options,1);
%% refine estimates excluding rejected components
%[A2,b2,C2] = update_spatial_components(Yr,C,f,[A,b],P,options);
%[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,P,options);
%% Plot refined ROI
%plot_contours(A2,Cn_max,options,1);
%% Do not refine estimates for soma:
refine=0;
if refine==0
A2=A;
b2=b;
C2=C;
P2=P;
f2=f;
S2=S;
YrA2=YrA;
end
%% Extract DF/F
[C_df,~] = extract_DF_F(Yr,A2,C2,P2,options);
% center for ROI 
center= com(A2,d1,d2);
% New df/f extract
alpha=0.05;
[expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA2, A2,C2, b2,f2,alpha); % extract DF/F values (optional)
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
%plot_components_GUI(Yr,A2,C2,b2,f2,options);
%% Save
tic;
save(fullfile([nam(1:end-6),session, '_Cdf']),'S','F_dff','C_df','expDffMedZeroed');
save(fullfile([nam(1:end-6),session,'_ROI']),'Coor' ,'json_file','center','A2','C2','b2','f2','Cn_max', 'options');
toc;

end
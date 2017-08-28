
%% Plot ROI
addpath(genpath('/Users/martial/Documents/GitHub/Phoenix'));
addpath(genpath('/Users/martial/Documents/GitHub/Canalysis'));

%% Open mat files with ROI

%For multiple files
listmat = dir('*.mat'); %will look for all the mat files
for i = 1:length(listmat)
ROI{i} = load((listmat(i).name));
end


%% Import tif image or use Cn
import_tif=0;
if import_tif==1
tif_file='/Users/martial/Desktop/PSAM_ROI/M991/STD_1_M991_CTRL._Tsub_mean.tif';
Cn_max = read_file(tif_file);
end
%% Options
% Plot only keep ROI:
plot_keep=1; % 0 or 1
plot_fr=10;
options.plot_df=1;
C_df=expDffMedZeroed';
%% plot
Y_plot =double( repmat(Cn_max,[1,1,plot_fr]) );
ROI=reshape(full(A2),length(Cn_max),length(Cn_max),size(A2,2));
All_ROI=sum(ROI,3);
figure;
imagesc(All_ROI)

if plot_keep==1 
plot_contours(A2(:,keep),Cn_max,options,1); 
plot_components_MAD(Y_plot,A2(:,keep),options,C_df);
elseif plot_keep==0
plot_contours(A2,Cn_max,options,1); 
plot_components_MAD(Y_plot,A2,options,C_df);
end


%% To try 
% To get the coordinates of dendrite ROIs
figure;
 imagesc(Cn_max);
      hold on;
  for r=1:size(ROI,3)
      for k = 1:length(B{r})
 %visboundaries(ROI(:,:,r))
 [B{r},L{r}] = bwboundaries(ROI(:,:,r));
   boundary = B{r}{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1)
coodonates{r}=B{r}{1};
end
end
  

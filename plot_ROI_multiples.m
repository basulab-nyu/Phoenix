
%% Plot ROI
addpath(genpath('/Users/martial/Documents/GitHub/Phoenix'));
addpath(genpath('/Users/martial/Documents/GitHub/Canalysis'));
%% Open mat files with ROI
listmat = dir('*.mat'); %will look for all the mat files
for i = 1:length(listmat)
ROI{i} = load((listmat(i).name));
end
%% Open mat files with cdf
listmat = dir('*.mat'); %will look for all the mat files
for i = 1:length(listmat)
cdf{i} = load((listmat(i).name));
C_df{i}=cdf{i}.expDffMedZeroed;
C_df{i}=(full(C_df{i})');
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
%% plot
FOV=length(ROI{1}.Cn_max);
if plot_keep==true; % 0 or 1
    if length(ROI)==3
keep_all=intersect(find(ROI{1}.keep==1), intersect(find(ROI{2}.keep==1), find(ROI{3}.keep==1)));
elseif length(ROI)==2
keep_all=intersect(find(ROI{1}.keep==1),find(ROI{2}.keep==1));
elseif length(ROI)==1
    keep_all=ROI.keep;
    end
elseif plot_keep==false;
    keep_all=1:size(C_df{1},2);
end
figure
for i=1:length(ROI)
Y_plot{i} =double( repmat(ROI{i}.Cn_max,[1,1,plot_fr]) );
A2_keep{i}=ROI{i}.A2(:,keep_all);
ROI_sh{i}=reshape(full(A2_keep{i}),FOV,FOV,size(A2_keep{i},2));
All_ROI{i}=sum(ROI_sh{i},3);
subplot(length(ROI),1,i)
imagesc(All_ROI{i})
end
for i=1:length(ROI)
    figure
plot_contours(A2_keep{i},ROI{i}.Cn_max,ROI{i}.options,1); 
legend(['ROI session ' num2str(i)]);
title(['session ' num2str(i)]);
end
for i=1:length(ROI)
Y_temp=Y_plot{i};
A2_temp=A2_keep{i};
options_temp=ROI{i}.options;
C_df_temp=C_df{i};
plot_components_MAD(Y_temp,A2_temp,options_temp,C_df_temp);
end



%% To try 
% To get the coordinates of dendrite ROIs
%set color of contour for sessions
color=[{'g'} {'r'} {'w'}];

for i=1:length(ROI)
Cn_all(:,:,i)=ROI{i}.Cn_max;
Cn_mean=mean(Cn_all,3);
leg{i}=['session ' num2str(i)];
end
figure
 imshow(Cn_mean);
      hold on;     
for i=1:length(ROI)
ROI_temp=ROI_sh{i};
for r=1:size(ROI_temp,3)
 %visboundaries(ROI(:,:,r))
[B{r},L{r}] = bwboundaries(ROI_temp(:,:,r));
for k = 1:length(B{r})
   boundary = B{r}{k};
   plot(boundary(:,2), boundary(:,1), color{i}, 'LineWidth', 1)
coodonates{r}=B{r}{1};

end
end
end 

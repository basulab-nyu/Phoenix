function plot_components_MAD(Y,A,options,C_df)

% GUI for plotting components that span large portions of the field of view (e.g., dendritic axonal)

memmaped = isobject(Y);
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
if ~isfield(options,'full_A') || isempty(options.full_A); full_A = defoptions.full_A; else full_A = options.full_A; end

if nnz(A)/numel(A) > 0.3; A = full(A); end
nA = full(sqrt(sum(A.^2))');
[K] = size(A,2);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy

nb=1;
Df=C_df;
fig = figure('Visible','off');
set(gcf,'Position',2*[300,300,960,480]);
set(gcf,'PaperPosition',2*[300,300,960,480]);

% Create a figure and axes
% ax = axes('Units','DF/F');
% Create slider
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',K+nb,'Value',1,'SliderStep',[1/(K+nb-1) 1],...
    'Position', [150 20 800 20],...
    'Callback', @surfzlim);

% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

% Make figure visble after adding all components
fig.Visible = 'on';
plot_component(1)

% This code uses dot notation to set properties.
% Dot notation runs in R2014b and later.
% For R2014a and earlier: set(f,'Visible','on');

    function surfzlim(source,callbackdata)
        i = source.Value;
        plot_component(round(i))
        % For R2014a and earlier:
        % i = get(source,'Value');               
    end

    function plot_component(i)
        subplot(2,2,[1,3]);
  %      if i <= K
            cla
            imagesc(reshape(A(:,i),d1,d2)); axis equal; axis tight; axis off; hold on;
  %          title(sprintf('Component %i ',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
  %      else
  %          cla
  %          imagesc(reshape(b(:,i-K),d1,d2)); axis equal; axis tight;
  %          title('Background component','fontsize',16,'fontweight','bold'); drawnow;
  %      end
        subplot(2,2,[2,4]);
   %     if i <= K
            plot(Df(:,i),'linewidth',2); hold all; plot(Df(:,i),'linewidth',2);
             title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
            drawnow;
            hold off;
   %     end 
    end
end
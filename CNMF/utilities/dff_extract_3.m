function [expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA2, A2,C2, b2,f2,alpha);

%non-denoised fluorescence for each component
% YrA2_nonL = A2'*Yr - (A2'*[A2,b2])*[C2;f2];
% A2_norm = 1./sum(A2.^2);
% A2_norm_rep = repmat(A2_norm',1,size(Yr, 2));
% 
% YrA2 = bsxfun(@times, YrA2_nonL, A2_norm_rep);

F = YrA2 + C2;

%background for non-denoised component
dfc = prctfilt(F,5,100,[],2);

%background (noise) fluorescence for each neuron
bf = prctfilt((bsxfun(@times, A2, 1./sum(A2.^2))'*b2)*f2,30,1000,300,2);

%dF/F background total (dfc + bf) = background + noise
dff = (F-dfc)./(dfc + bf);

%exponential moving average to smooth
expDff = []; 

for idx=1:size(dff,1)
    expDff(idx,:) = filter(alpha, [1 alpha-1], dff(idx,:));
end

%substract median from each trace to baseline around zero
for idx=1:size(expDff,1)
    expDffMedZeroed(idx,:) = expDff(idx,:) - median(expDff(idx,:));
end

end
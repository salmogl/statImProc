clc;
clear all;
close all;

% This script applies image denoising with 3 different operators: 
% A             : regular diffusion operator (rows normalization, NLM)
% A_            : diffusion operator (columns and rows normalization)
% 2A_ - A_^2    : diffusion operator suggested in (1) (columns and rows normalization)
%
% Relevant Papers:
% (1) "DIFFUSION INTERPRETATION OF  NON-LOCAL NEIGHBORHOOD FILTERS FOR SIGNAL
%     DENOISING [AMIT SINGER, YOEL SHKOLNISKY, BOAZ NADLER]"
% (2) "Demystifying Symmetric Smoothing Filters [Stanley H. Chan, Todd
%     Zickler, Yue M. Lu]
%==========================================================================
% AUTHOR        Almog Lahav
% INSTITUTION   Technion
% DATE          23th August 2016
%
%
% SCRIPT PARAMETERS  see inside the script
%==========================================================================

Im = double(imread('barbara.png'));

Im = Im(80:128,382+15:430+15); % 49X49 pxl
% Im = Im(29:128,318+28:445); % 100X100 pxl
% Im = Im(1:128,318:445); % 128X128 pxl


Im = Im - mean(Im(:));
Im = Im/std(Im(:));
patchSize = 5;

margin = round(patchSize/2);
idxImNoMarginRows = margin:size(Im,1)-margin+1;
idxImNoMarginCols = margin:size(Im,2)-margin+1;
ImTrunc = Im(idxImNoMarginRows,idxImNoMarginCols);

stdNoise = 60/255;

params.normCols = 2;          % 1 = columns and rows normalization; 2 = rows normalization only
params.metric = 'euc';        % metric for the kernel
params.knn = 30;              % number of nearest neighboors (nn) to compute the kernel bandwidth
params.eps = 0.5;             % fraction of the median of the nn
params.thresh = 1e-8;         % under params.thresh the kernel is 0
params.freqfilt = false;      % filtering using SVD


% For optimal bandwidth run the script comparePsnr.m
eps = [0.425 0.425 0.65]; % optimal bandwidh for the sub image 128x128.  
pSnrPatch = zeros(1,3);
allImFiltered = zeros([size(ImTrunc),3]);

ImNoisy = Im + stdNoise * randn(size(Im));
[patches,ImAllPatches] = getPatches(ImNoisy,patchSize);
ImNoisyTrunc = ImNoisy(idxImNoMarginRows,idxImNoMarginCols);
       
for j = 1:3
    
    params.eps = eps(j);
    % alternating between normalized and un-normalized
    params.normCols = mod(j-1,2)+1;

        [ImFilteredPatch,eigVectPatch,eigValsPatch,APatch] = diffusionFilter(patches,ImAllPatches,params);
    if j == 2 
        ImFilteredPatch = reshape((2*APatch - APatch^2)*ImAllPatches,size(ImTrunc));
    end
    pSnrPatch(j) = psnr(ImFilteredPatch,ImTrunc);
    pSnrNosiy = psnr(ImNoisyTrunc,ImTrunc);
    allImFiltered(:,:,j) = ImFilteredPatch;

end



figure;
subplot(221)
imshow(ImNoisyTrunc,[])
title('(a) Noisy Image','Interpreter','latex','Fontsize',15)
subplot(222)
imshow(allImFiltered(:,:,1),[])
title(strcat('(b) $\tilde{A}$, PSNR = ',num2str(pSnrPatch(1)))...
    ,'Interpreter','latex','Fontsize',15)
subplot(223)
imshow(allImFiltered(:,:,2),[])
title(strcat('(c) $A$, PSNR = ',num2str(pSnrPatch(2)))...
    ,'Interpreter','latex','Fontsize',15)
subplot(224)
imshow(allImFiltered(:,:,3),[])
title(strcat('(d) $2\tilde{A}-\tilde{A}^2$, PSNR = ',num2str(pSnrPatch(3)))...
    ,'Interpreter','latex','Fontsize',15)

figure;
subplot(221)
imshow(ImNoisyTrunc-allImFiltered(:,:,2),[])
title('(a) $I_n-I_A$','Interpreter','latex','Fontsize',15)
subplot(222)
imshow(ImNoisyTrunc-allImFiltered(:,:,3),[])
title('(b)  $I_n-I_{\tilde{A_2}}$','Interpreter','latex','Fontsize',15)
subplot(223)
imshow(ImTrunc-allImFiltered(:,:,2),[])
title('(c)  $I-I_A$','Interpreter','latex','Fontsize',15)
subplot(224)
imshow(ImTrunc-allImFiltered(:,:,3),[])
title('(d)  $I-I_{\tilde{A_2}}$','Interpreter','latex','Fontsize',15)


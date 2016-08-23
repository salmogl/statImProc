clc;
clear all;
close all;

% This script compares the psnr of image denoising with 3 different operators: 
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

maxIter = 1;
eps = linspace(0.4,1.2,15);% linspace(0.005,1.2,5); %linspace(0.005,1.2,5); %0.667

pSnrPatch = zeros(maxIter,length(eps),3);


for iter = 1:maxIter

    ImNoisy = Im + stdNoise * randn(size(Im));
    [patches,ImAllPatches] = getPatches(ImNoisy,patchSize);
    

    ImNoisyTrunc = ImNoisy(idxImNoMarginRows,idxImNoMarginCols);

    for iEps = 1:length(eps)
        
        params.eps = eps(iEps);
        
        for j = 1:3
            
            % alternating between normalized and un-normalized
            params.normCols = mod(j-1,2)+1;

                [ImFilteredPatch,eigVectPatch,eigValsPatch,APatch] = diffusionFilter(patches,ImAllPatches,params);
            if j == 3 
                ImFilteredPatch = reshape((2*APatch - APatch^2)*ImAllPatches,size(ImTrunc));
            end
            pSnrPatch(iter,iEps,j) = psnr(ImFilteredPatch,ImTrunc);
            pSnrNosiy = psnr(ImNoisyTrunc,ImTrunc);
            
            
        end

    end
end

meanPsnrPatch = squeeze(mean(pSnrPatch,1));

figure;
plot(eps,meanPsnrPatch')
legend('Columns & Rows Normalization','Rows Normalization','2A-A^2')


figure;
imshow(ImTrunc,[])

figure;
imshow(ImNoisyTrunc,[])



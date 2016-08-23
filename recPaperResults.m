clc;
clear all;
close all;

% This script recovers results from the paper:
% "DIFFUSION INTERPRETATION OF  NON-LOCAL NEIGHBORHOOD FILTERS FOR SIGNAL
% DENOISING [AMIT SINGER, YOEL SHKOLNISKY, BOAZ NADLER]"
%==========================================================================
% AUTHOR        Almog Lahav
% INSTITUTION   Technion
% DATE          23th August 2016
%
%
% SCRIPT PARAMETERS  see inside the script
%==========================================================================

Im = double(imread('barbara.png'));

% get sub-images
Im = Im(80:128,382+15:430+15); % 49X49   pxl
% Im = Im(29:128,318+28:445);  % 100X100 pxl
% Im = Im(1:128,318:445);      % 128X128 pxl


Im = Im - mean(Im(:));
Im = Im/std(Im(:));
stdNoise = 60/255;
ImNoisy = Im + stdNoise * randn(size(Im));

patchSize = 5;
[patches,ImAllPatches] = getPatches(ImNoisy,patchSize);

% remove margins 
margin = round(patchSize/2);
idxImNoMarginRows = margin:size(Im,1)-margin+1;
idxImNoMarginCols = margin:size(Im,2)-margin+1;
ImTrunc = Im(idxImNoMarginRows,idxImNoMarginCols);
ImNoisyTrunc = ImNoisy(idxImNoMarginRows,idxImNoMarginCols);

coordinates = getCoordinates(ImNoisyTrunc);

params.normCols = 2;          % 1 = columns and rows normalization; 2 = rows normalization only
params.metric = 'euc';        % metric for the kernel
params.knn = 30;              % number of nearest neighboors (nn) to compute the kernel bandwidth
params.eps = 0.5;             % fraction of the median of the nn
params.thresh = 1e-8;         % under params.thresh the kernel is 0
params.freqfilt = true;       % filtering using SVD


% filter the image using kernel of patches and coordinates
% [ImFilteredPatch,eigVectPatch,eigValsPatch,A] = diffusionFilter([patches;coordinates],ImAllPatches,params);

% filter the image using kernel of patches
[ImFilteredPatch,eigVectPatch,eigValsPatch,APatch] = diffusionFilter(patches,ImAllPatches,params);

% filter the image using kernel of single pixles
params.eps = 0.3;
[ImFilteredCoor,eigVectCoor,eigValsCoor,ACoor] = diffusionFilter(coordinates,ImNoisyTrunc(:),params);


% compute errors
errPatch = norm(ImTrunc-ImFilteredPatch,'fro');
errCoor = norm(ImTrunc-ImFilteredCoor,'fro');
pSnrCoor = psnr(ImFilteredCoor,ImTrunc);
pSnrPatch = psnr(ImFilteredPatch,ImTrunc);
pSnrNosiy = psnr(ImNoisyTrunc,ImTrunc);


% plot images
figure;
subplot(221)
imshow(ImTrunc,[])
title('Barbara (Original)')

subplot(222)
imshow(ImNoisyTrunc,[])
title('Barbara (Noisy)')

subplot(223)
imshow(ImFilteredPatch ,[])
title('Barbara (Denoised - Patch Diffusion)')

subplot(224)
imshow(ImFilteredCoor,[])
title('Barbara (Denoised - Coordinates Diffusion)')

figure;
subplot(121)
plot(eigValsCoor)
title('Eignvalues Of The Coordinates Diffusion Kernel')

subplot(122)
plot(eigValsPatch)
title('Eignvalues Of The Patches Diffusion Kernel')

% plot the embedding of the data in the diffusion space. Color points by
% their transition probability from them to the point: pointIndex
 if (params.freqfilt)
    
    pointIndex = 309;
    nRadius = 50;
    prob = linspace(0,1,nRadius).^15;
    colors = flipud(jet(nRadius));
    C  = zeros(length(eigVectCoor),3);
    nighbors = APatch(pointIndex,:);

    for i = 1:nRadius-1
        M = length(find(nighbors < prob(i+1) & nighbors >= prob(i)));
        C(nighbors < prob(i+1) & nighbors >= prob(i),:) = repmat(colors(i,:),M,1);
    end


    C(pointIndex,:) = zeros(1,3);
    N = length(eigVectCoor);
    embedding = eigVectPatch(:,2:4)*diag(eigValsPatch(2:4));

    figure;
    scatter3(embedding(:,1),embedding(:,2),embedding(:,3),[],C,'filled');
    xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3')
    text(embedding(pointIndex,1),embedding(pointIndex,2),embedding(pointIndex,3)+0.03,'\downarrow','fontsize',40)
    title('Embedding Using Patches Diffusion Kernel')
    caxis([0 1])
    h = colorbar();
    set(h,'yDir','reverse')
    colormap jet

 end



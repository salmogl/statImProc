function [ImFiltered,eigVect,eigVals,A] = diffusionFilter(patches,ImAllPatches,params)
% diffusionFilter  filters the image ImAllPatches using a diffusion
% operator (=NLM)
%==========================================================================
% AUTHOR        Almog Lahav
% INSTITUTION   Technion
% DATE          23th Auguat 2016
%
%
% INPUTS
%   patches                   - all the patches of the image. Full overlapping. Each
%                               column is patch
%   ImAllPatches              - columns stack of all pixels not including margins.
%   params.normCols           - 1 = columns and rows normalization; 2 = rows normalization only
%   params.metric             - metric for the kernel
%   params.knn                - number of nearest neighboors (nn) to compute the kernel bandwidth
%   params.eps                - fraction of the median of the nn
%   params.thresh             - under params.thresh the kernel is 0
%   params.freqfilt           - filtering using SVD if true
%
% OUTPUTS
%   ImFiltered   - The filtered image
%   eigVect      - Left eigenvectors of the diffusion operator. 0 if
%                  params.freqfilt == false
%   eigVals      - Eigenvalues of the diffusion operator. 0 if
%                  params.freqfilt == false
%   A            - The diffusion operator
%==========================================================================

 W = CalcInitAff( patches, params );

switch params.normCols
    
    case 1
        
        Dc = diag(sum(W,1));
        Dr = diag(sum(W*Dc^-1,2));
        A = (Dr^-1)*W*(Dc^-1);
        
    case 2
                
        D = diag(sum(W,2));
        A = D^-1*W;
        
    case 3
        
        eucDist = squareform(pdist(patches'));
        [~, nn_dist] = knnsearch(patches', patches', 'k', params.knn);
        sigma = params.eps * median(nn_dist(:));
        W = exp(-eucDist.^2./(2*sigma))./sqrt(2*pi*sigma);
        
        P0 = repmat(mean(W,1),size(W,1),1);
        G = geomean(mean(W,1));
        epsj = (G./P0);
        W_ = exp(-eucDist.^2./(2*sigma*epsj))./sqrt(2*pi*sigma*epsj);
        D = diag(sum(W_,2));
        A = D^-1*W;
        
    case 4
                
        D = diag(sum(W,1));
        W0 = D*W*D^-1;
        Dr = diag(sum(W0,2));
        A = Dr^-1*W0;
        
    case 5
                
        Dr = diag(sum(W,2));
        Dc = diag(sum(W,1));
        W = (Dr^-0.5)*W*(Dc-0.5);
        D = diag(sum(W,2));
        A = D^-1*W;
    
           
end
 
N = sqrt(length(ImAllPatches)); 

if (params.freqfilt)
    [eigVect,S,V] = svd(A);
    % eigVect =eigVect(:,1:4);
    % S = S(1:4,1:4);
    % V =V(:,1:4);
    eigVals = sum(S);
    ImFiltered = eigVect*S*V'*ImAllPatches;
    ImFiltered = reshape(ImFiltered,N,N);
else
    eigVals = zeros(1,N^2);
    eigVect = zeros(N^2,N^2);
    ImFiltered = reshape(A*ImAllPatches,N,N);
end
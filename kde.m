function [p,yi] = kde(y,args)
% kde  computes the density estimation p in the points yi
%==========================================================================
% AUTHOR        Almog Lahav
% INSTITUTION   Technion
% DATE          23th Auguat 2016
%
%
% INPUTS
%   y            - data. Columns ar samples, rows are features.
%   args.norm    - methode for density estimation:
%                  1 = fixed kernel. (equivalent to rows normalization)
%                  2 = normalized kernel (equivalent to columns and rows normalization)
%                  3 = addaptive kernel (equivalent to Addaptive KDE)
% OUTPUTS
%   p            - values of the estimated density at the points yi
%   yi           - the points in which the density is calculated
%==========================================================================

eucDist = squareform(pdist(y'));
[IDX, nn_dist] = knnsearch(y', y', 'k', args.knn);
sigma = args.eps*median(nn_dist(:));
W = exp(-eucDist.^2/(2*sigma))/sqrt(2*pi*sigma);

switch args.norm
    
    case 'rows'
    
        D = mean(W,2)';
    
    case 'columns'

        Dc = diag(mean(W,1));
        D = mean(W*Dc^-1,2)'.*mean(W,1);
    
    case 'adaptive'
        
        P0 = repmat(mean(W,1),length(y),1);
        G = geomean(mean(W,1));
        epsj = (G./P0);
        W = exp(-eucDist.^2./(2*sigma*epsj))./sqrt(2*pi*sigma*epsj);
        D = mean(W,2)';
end

[yi,index] = sort(y);
p = D(index);
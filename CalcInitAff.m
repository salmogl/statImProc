function [ aff_mat ] = CalcInitAff( data, params )
% CalcInitAff  Calculates the affinity between the columns of data
%==========================================================================
% INPUTS
%   data                      - M-by-N matrix whose columns are N points in R^M 
%   params.metric             - metric for the kernel
%   params.knn                - number of nearest neighboors (nn) to compute the kernel bandwidth
%   params.eps                - fraction of the median of the knn
%   params.thresh             - under params.thresh the kernel is 0
%
% OUTPUTS
%   aff_mat                   - Affinity matrix
%==========================================================================
switch params.metric
    case 'cosine_similarity'
    inner_products = data.'*data;
    [ij_norm, ji_norm] = meshgrid(diag(inner_products));
    aff_mat = inner_products./sqrt(ij_norm.*ji_norm);
    aff_mat(aff_mat < params.thresh) = 0;
    case 'euc'
    data = data.';
%     C = diag([(10*params.stdNoise)^2*(ones(1,size(data,2)-2)) params.eps^2 params.eps^2]);
%     euc_dist = squareform(pdist(data,'mahalanobis',C));
    euc_dist = squareform(pdist(data));

    [IDX, nn_dist] = knnsearch(data, data, 'k', params.knn);
    sigma = params.eps * median(nn_dist(:));
    aff_mat = (2*pi*sigma)^-0.5*exp(-(euc_dist.^2)/(2*sigma)); 
    aff_mat(aff_mat < params.thresh) = 0;
           
        
end

% figure, imagesc(aff_mat), colormap jet, colorbar, axis on, hold on
% if params.on_rows,
%     title('Initial Row Affinity'), hold off
% else
%     title('Initial Col Affinity'), hold off
% end

end


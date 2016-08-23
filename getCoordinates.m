function [coordinates] = getCoordinates(Im)
% getCoordinates return the spatial coordinates of all pixels of the image
% Im
%==========================================================================
% INPUTS
%   Im                        - Image 

% OUTPUTS
%   coordinates               - Spatial coordinates. Each column in 
%                               coordinates is a pair of spatial coordinates
%==========================================================================


[M,N] = size(Im);

coordinates = zeros(2,M*N);
patchIndex = 1;

for j = 1:N
    for i = 1:M
        
        coordinates(:,patchIndex) = [i;j];
        patchIndex = patchIndex + 1;
        
    end
end


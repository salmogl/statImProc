function [patches,ImAllPatches] = getPatches(Im,patchSize)

% getPatches returns all the patches of the size patchSize X patchSize in
% Im. Full overlaping.
%==========================================================================
% INPUTS
%   Im                        - Image 
%   patchSize                 - Patch size
%
% OUTPUTS
%   patches                   - Each column is a patch 
%   ImAllPatches              - All the pixels of Im that has patches (does not
%                               include the margin) in column stack;
%==========================================================================
 
patches = zeros(patchSize^2,(size(Im,1)-patchSize+1)*(size(Im,2)-patchSize+1));
patchIndex = 1;

for j = 1:size(Im,2)-patchSize+1
    for i = 1:size(Im,1)-patchSize+1 
        
        curPatch = Im(i:i+patchSize-1,j:j+patchSize-1);
        patches(:,patchIndex) = curPatch(:);
        patchIndex = patchIndex + 1;
        
    end
end

ImAllPatches = Im(round(patchSize/2):end-round(patchSize/2)+1,...
    round(patchSize/2):end-round(patchSize/2)+1);
ImAllPatches = ImAllPatches(:);
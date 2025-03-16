function newMaskImg = combineMask(mask1, mask2, newMask)
%newMask = combineMask(mask1, mask2, newMask)
% mask1, mask2, newMask: bmp mask names

if(nargin==2)
    newMask = 'new_SAXS_mask.bmp';
elseif (nargin<2)
    fprintf('Not enough parameters! Usage: combineMask(mask1, mask2, newMask)');
end

% mask value: 0 belongs mask; 1 belongs image
m1 = logical(imread(mask1));
m2 = logical(imread(mask2));
newMaskImg = logical(m1.*m2);
imwrite(newMaskImg, newMask);
end


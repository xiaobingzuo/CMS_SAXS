function genNewPixelMask(dkFileTemp, num, maskFile)
%genNewPixelMask(dkFileTemp, num)
%Usage: 
% (1)in spec: takeshot dark 0.1 10 2
% (2)Assuming exact dark files are: SAXS\Sdark_00032_00001.h5 - _00010.h5  
%    in Matlab: genNewPixelMask('SAXS\Sdark_00032', 10)
% (3)'new_pixel_mask.bmp' will be generated for mask use. Dead pixel# after
%     factory calibration will be printed out.
%
%Note: Inactive pixel number when leaving factory: 664597 ( 664248 in module gaps, and 349
%dead pixels
if (nargin == 2)
    maskFile = 'new_pixel_mask.bmp';
elseif (nargin < 2)
    disp('Not enough parameters, Usage: genNewPixelMask(darkFileTemp, #ofDarkImage, maskFile');
    return;
end

mask9 = [];
for kk = 1:num
    cfn = sprintf('%s_%05d.h5', dkFileTemp, kk);
    mask0 = readHdf5(cfn);
    cmask = logical(mask0.data');
    if kk==1
        mask9 = cmask;
    else
        mask9 = mask9 .* cmask;
    end
end
imwrite(logical(flipud(1-mask9)), maskFile);
figure(15); imagesc(1-mask9);
fprintf('New dead pixel#: %d.\n', sum(sum(mask9>0)) - 664597);
end


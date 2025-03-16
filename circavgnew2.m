function data = circavgnew2(imag, mask, qCorrMap, qRMap, qArray, offset, limits)
% data = circavgnew2(image, mask, qCorrMap, qRMap, qArray, offset, Limits)
% Inputs: 
%    image:    image in matrix to be process
%    mask:     pre-made mask file in matrix
%    qCorrMap: pre-calculated solid angle correction map for pixels on detector
%    qArray:   pre-defined q-value array
%    qRMap:    q-index map of pixels on the detector for the given qArray 
%    offset:   intensity offset value 
%    Limits:   the lowest and highest count values will be considered as valid counts
% Output:
%    data: 1-d scattering profile converted from the 2D image. 'data' contains five columns:
%          q, Iq, Iq_std, tot_corted_count_at_q, tot_valid_pixel_num, Iq_SA-uncorrected

% Written by Xiaobing Zuo, 03/2011

data=[];

if isempty(mask)
    mask = ones(size(imag));
end

imgSize=size(imag);
mskSize=size(mask);
qcMapSize=size(qCorrMap);
qrMapSize=size(qRMap);
rows=imgSize(1);
cols=imgSize(2);
if imgSize==mskSize & imgSize==qcMapSize & imgSize==qrMapSize
     %data = circularavgnew(double(imag), double(mask), qCorrMap, qRMap, qArray, offset, limits);
     %data = circularavgnew2(double(imag), double(mask), qCorrMap, qRMap, qArray, offset, limits);

     data = circularavgnew2b(double(imag), double(mask), qCorrMap, qRMap, qArray, offset, limits);
     % data cols: 1_q   2_iq  3_sqrt(iq/N)   4_mean(iq_i^2) 5_pixel#  6_mean(corr_i)  
else
    ['image: (', num2str(imgSize), '); mask: (',  num2str(mskSize), '); qcorrMap: (', num2str(qcMapSize), '); qRMap: (', num2str(qrMapSize), ').']
end    


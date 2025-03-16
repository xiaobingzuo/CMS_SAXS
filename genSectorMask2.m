function genSectorMask2(voxels,col,row)
%[qCalData, userData, scatt, calibrant]= getData;
tic
imgx = [];
imgy = [];
xv = voxels(:,1); yv =voxels(:,2);
min(xv)
max(xv)
min(yv)
max(yv)
%col = scatt.col; row = scatt.row;
for kk=1:col
    imgx = [imgx ones(1,row)*kk];
    imgy = [imgy 1:row];
end
%imgx(:).'
toc
tic
in = inpolygon(imgx,imgy,xv',yv');
mask2 = reshape(in,[row, col]);
toc
figure; imagesc(mask2);
figure; imagesc(1-mask2);
% if scatt.sectorMaskon
%     if scatt.sectorMaskExclude
%         mask2 = [];
%     else
%         mask2 = [];
%     end
% else
%     mask2 = [];
% end
%setData(qCalData, userData, scatt, calibrant);
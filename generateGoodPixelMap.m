function mask_inUse = generateGoodPixelMap(fn, threshold, transp, folder)
%This function generate goodpix_mask_c.bmp
%how to use: take a series of images without X-ray
% e.g.; takeshot blank 1 10 2; to generate 10 '*blank_00*'
%  mask_inUse = generateGoodPixelMap(fn, threshold,  transp, folder)
%  fn: file name to process, in above case use 'blank'
%  threshold: without x-ray, the count should be 0, threshold can be 1 to
%        remove cosmic ray counts
%  transp: 1 if need to transpose image first, Eiger data at 12-ID-B needs
%           to transpose; Pilatus image does not need, should be 0. 
%  folder: the folder where those image located
%
% example commands:
% For 12-ID-B Eiger:
% In spec: takeshot blank 1 10 2; % to generate 10 'blank_00001_*.h5'
% generateGoodPixelMap('blank_00001', 2, 1)

if nargin == 1    
    threshold =2;
    transp=1;
    folder = pwd;
elseif  nargin == 2
    transp=1;
    folder = pwd;      
elseif  nargin == 3
    folder = pwd;    
elseif (nargin <1 || nargin >3) 
    fprintf('input arguments wrong!');
end

fn0=['*' fn '*'];
fn1 = fullfile(folder, fn0);

mask=[];

fn2 = dir(fn1);

if ~isempty(fn2)
    h5=readHdf5(fullfile( fn2(1).folder, fn2(1).name), '12-ID-B');
    if transp
        mask=(h5.data>= threshold)';   
    else
        mask=(h5.data>= threshold); 
    end
    figure(2); imagesc(mask);
    [1, sum(sum(mask))]
    for jj=1:numel(fn2)
        h5=readHdf5(fullfile( fn2(jj).folder, fn2(jj).name), '12-ID-B');
        if transp
            cMask = (h5.data>= threshold)';
        else
            cMask = (h5.data>= threshold);
        end
        mask = mask & cMask;
        [jj, sum(sum(cMask)), sum(sum(mask))]
    end
end
figure(3); imagesc(mask);
sum(sum(mask))
figure(4); imagesc(flipud(mask));
mask_inUse = flipud(~mask);
imwrite(mask_inUse, fullfile(folder, 'goodpix_mask_c.bmp'));
end
function  readoutBNLraw1D(h5FileName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

global fileIndex;

fnInfo = h5info(h5FileName);
nSamp = numel(fnInfo.Groups);
cfolder = pwd;
for kk=1:nSamp
    cfn0 = fnInfo.Groups(kk).Name;
    cfn =cfn0(2:end);
    qgrid = h5readatt(h5FileName,'/','qgrid');
    %dat1=h5read(fn, '/css7/processed/_SAXS'); 
    %dat2=h5read(fn, '/css7/processed/_WAXS2'); 
    dat3=h5read(h5FileName, [cfn0 '/processed/averaged']);
    dat4=h5read(h5FileName, [cfn0 '/processed/merged']);
    
    processedFolder =fullfile(cfolder, 'Processed');
    if ~exist(processedFolder, 'dir')
       mkdir(processedFolder)
    end

    AveragedFolder =fullfile(cfolder, 'Averaged');
    if ~exist(AveragedFolder, 'dir')
       mkdir(AveragedFolder)
    end
    

    %fnavg=[cfn '_a.dat'];
    fnavg1 = sprintf('%s_%05d_avg.dat', cfn, fileIndex);
    %fnavg =fullfile(cfolder, 'Processed', fnavg1);
    fnavg =fullfile(processedFolder, fnavg1);
    dlmwrite(fnavg, [qgrid dat3], 'delimiter', ' ', 'precision', '%.4E');

    try
        dat5=h5read(h5FileName, [cfn0 '/processed/subtracted']);
        %fnsub=[cfn '_s.dat'];
        fnsub1 = sprintf('%s_%05d_sub.dat', cfn, fileIndex);
        fnsub= fullfile(processedFolder, fnsub1);
        %dlmwrite(fullfile(resultFolder, fnsub), [qgrid dat5], 'delimiter', ' ', 'precision', '%.4E');
        dlmwrite(fnsub, [qgrid dat5], 'delimiter', ' ', 'precision', '%.4E');
    catch
        errstr = sprintf('File: %s, has no subtracted file!', cfn);
        disp(errstr);
    end

    [~,~,nImg] = size(dat4);
    for kk = 1:nImg
        fn1a = sprintf('%s_%05d_%05d.dat', cfn, fileIndex, kk);
        fn1 = fullfile(AveragedFolder, fn1a);
        dlmwrite(fn1, [qgrid dat4(:,:,kk)], 'delimiter', ' ', 'precision', '%.4E');
    end
    fileIndex = fileIndex + 1;
end


end
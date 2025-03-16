function  readoutBNL1D(h5FileName,resultFolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

%global fileIndex;

fnInfo = h5info(h5FileName);
nSamp = numel(fnInfo.Groups);

for kk=1:nSamp
    cfn0 = fnInfo.Groups(kk).Name;
    cfn =cfn0(2:end);
    qgrid = h5readatt(h5FileName,'/','qgrid');
    %dat1=h5read(fn, '/css7/processed/_SAXS'); 
    %dat2=h5read(fn, '/css7/processed/_WAXS2'); 
    dat3=h5read(h5FileName, [cfn0 '/processed/averaged']);
    dat4=h5read(h5FileName, [cfn0 '/processed/merged']);

    fnavg=[cfn '_a.dat'];
    dlmwrite(fullfile(resultFolder, fnavg), [qgrid dat3], 'delimiter', ' ', 'precision', '%.4E');

    try
        dat5=h5read(h5FileName, [cfn0 '/processed/subtracted']);
        fnsub=[cfn '_s.dat'];
        dlmwrite(fullfile(resultFolder, fnsub), [qgrid dat5], 'delimiter', ' ', 'precision', '%.4E');
    catch
        errstr = sprintf('File: %s, has no subtracted file!', cfn);
        disp(errstr);
    end

    %[~,~,nImg] = size(dat4);

end


end
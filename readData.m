function mat=readData(fName, nskip)
%  mat=readData(fName, nskip)
%  fName: name of file to read
%  nskip: n-lines to skip
%  Return the data into matrix mat.
%  
% This function returns a matrix of the data file while ignores non-data portion in the beginning of the file.

% by Xiaobing Zuo, Argonne National Laboratory
kk=1;
if nargin==1
    nskip = 0;
end

fid=fopen(fName);
if fid==-1
    disp(sprintf('%s was not able to open!',fName));
    mat = [];
    return ;
else
    mat=[];
    %nLine = 1;
    while feof(fid) ~= 1
        cLine = fgetl(fid);
        cdat = str2num(cLine);
        [nx ny] = size(cdat);
        if nx == 1 & kk > nskip
            mat = [mat; cdat ];
        end
        kk=kk+1;
    end    
    fclose(fid);
end
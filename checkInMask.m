function found = checkInMask(mt,mask)
% found = checkInMask(mt,mask)
% check if any point/item of mx in mask
% mask value == 0; ==1 data region
found = 0;
if ~isempty(mask)   % if empty, found =0
    for kk=1:length(mt)
        %[mt(kk,1), mt(kk,2)]
        if mask(mt(kk,2), mt(kk,1)) == 0
            found = 1;
            break;
        end
    end
end
end
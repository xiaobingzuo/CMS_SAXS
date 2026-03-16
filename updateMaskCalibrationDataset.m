function updateMaskCalibrationDataset(maskSets, calibSets, destFolder)
% updateMaskCalibrationDataset
% Manage mask and calibration datasets and record metadata in JSON.
%
% INPUTS
% maskSets  : struct array with fields
%             .maskFile
%             .motherFile
%             .type   ('SAXS' or 'WAXS')
%
% calibSets : struct array with fields
%             .imageFile
%             .setupFile
%             .type ('SAXS' or 'WAXS')
%
% destFolder : destination folder for storing datasets
%
% JSON file written: Mask_Calibration_dataset.json

MAX_Sets = 4;

jsonFile = fullfile(destFolder,'Mask_Calibration_dataset.json');

%% Create folder if needed
if ~exist(destFolder,'dir')
    mkdir(destFolder);
end

%% Load existing JSON dataSets
if exist(jsonFile,'file')
    fid = fopen(jsonFile);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    dataSets = jsondecode(str);
else
    dataSets.maskSets = struct([]);
    dataSets.calibrationSets = struct([]);
end

%% -------------------------
% PROCESS MASK DATASETS
%% -------------------------

for i = 1:length(maskSets)

    maskName = getFileName(maskSets(i).maskFileName);
    motherName = getFileName(maskSets(i).motherFileName);

    % check if identical mask already exists
    existingIdx = [];
    for j = 1:length(dataSets.maskSets)
        if strcmp(dataSets.maskSets(j).maskFileName,maskName)
            existingIdx = j;
            break
        end
    end

    % limit to 4 sets
    if isempty(existingIdx) && length(dataSets.maskSets) >= MAX_Sets
        error('Maximum of 4 mask datasets allowed.');
    end

    % copy files
    copyfile(fullfile(maskSets(i).maskOriginalLocation, maskSets(i).maskFileName), fullfile(destFolder, maskName));
    copyfile(fullfile(maskSets(i).motherOriginalLocation, maskSets(i).motherFileName), fullfile(destFolder, motherName));

    newEntry.maskFileName = maskName;
    newEntry.motherFileName = motherName;
    if isfield(maskSets(i), 'maskOriginalLocation')
        newEntry.maskOriginalLocation = maskSets(i).maskOriginalLocation;
    else
        newEntry.maskOriginalLocation = maskSets(i).maskFileName;
    end
    
    if isfield(maskSets(i), 'motherOriginalLocation')
        newEntry.motherOriginalLocation = maskSets(i).motherOriginalLocation;
    else
        newEntry.motherOriginalLocation = maskSets(i).motherFileName;
    end
    
    %newEntry.motherOriginalLocation = maskSets(i).motherFile;
    newEntry.type = maskSets(i).type;
    newEntry.timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');

    % overwrite or append
    if isempty(dataSets.maskSets)
        dataSets.maskSets = newEntry;
    else
        if isempty(existingIdx)
            dataSets.maskSets(end+1) = newEntry;
        else
            dataSets.maskSets(existingIdx) = newEntry;
        end
    end

end


%% -------------------------
% PROCESS CALIBRATION DATASETS
%% -------------------------

for i = 1:length(calibSets)

    setupName = getFileName(calibSets(i).setupFileName);
    imgName = getFileName(calibSets(i).imageFileName);

    % check uniqueness of setup file
    existsFlag = false;
    for j = 1:length(dataSets.calibrationSets)
        if strcmp(dataSets.calibrationSets(j).setupFileName,setupName)
            existsFlag = true;
            idx = j;
            break
        end
    end

    % limit to 4 --> max_sets
    if ~existsFlag && length(dataSets.calibrationSets) >=  MAX_Sets
        error('Maximum of 4 calibration datasets allowed.');
    end

    % copy files
    copyfile(fullfile(calibSets(i).imageOriginalLocation, calibSets(i).imageFileName), fullfile(destFolder,imgName));
    copyfile(fullfile(calibSets(i).setupOriginalLocation, calibSets(i).setupFileName), fullfile(destFolder,setupName));

    newCalibEntry.imageFileName = imgName;
    newCalibEntry.setupFileName = setupName;
    if isfield(calibSets(i), 'imageOriginalLocation')
        newCalibEntry.imageOriginalLocation = calibSets(i).imageOriginalLocation;
    else
        newCalibEntry.imageOriginalLocation = calibSets(i).imageFileName;
    end

    if isfield(calibSets(i), 'setupOriginalLocation' )
        newCalibEntry.setupOriginalLocation = calibSets(i).setupOriginalLocation;
    else
        newCalibEntry.setupOriginalLocation = calibSets(i).setupFileName;
    end
    newCalibEntry.type = calibSets(i).type;    
    newCalibEntry.calibrantName = calibSets(i).calibrantName;
    newCalibEntry.calibrantQValues = calibSets(i).calibrantQValues;
    newCalibEntry.XrayEnergy = calibSets(i).XrayEnergy;
    newCalibEntry.timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');

    % update dataSets
    if isempty(dataSets.calibrationSets)
        dataSets.calibrationSets = newCalibEntry;
    else
        if existsFlag
            dataSets.calibrationSets(idx) = newCalibEntry;
        else
            dataSets.calibrationSets(end+1) = newCalibEntry;
        end
    end

end


%% -------------------------
% SAVE JSON
%% -------------------------

jsonText = jsonencode(dataSets,'PrettyPrint',true);

fid = fopen(jsonFile,'w');
fwrite(fid,jsonText,'char');
fclose(fid);

disp('Mask/Calibration dataset successfully updated.');

end


%% Helper function
function name = getFileName(pathstr)
[~,name,ext] = fileparts(pathstr);
name = [name ext];
end

%% Example Usage
% maskSets(1).maskFile = 'C:\data\mask1.tif';
% maskSets(1).motherFile = 'C:\data\sample1.h5';
% maskSets(1).type = 'SAXS';
% 
% maskSets(2).maskFile = 'C:\data\mask2.tif';
% maskSets(2).motherFile = 'C:\data\sample2.tif';
% maskSets(2).type = 'WAXS';
% 
% 
% calibSets(1).imageFile = 'C:\data\calib_img1.tif';
% calibSets(1).setupFile = 'C:\data\calib_setup1.txt';
% calibSets(1).type = 'SAXS';
% 
% 
% updateMaskCalibrationDataset(maskSets, calibSets, 'Mask_Calibration_TrainingSets');
% Example JSON Output
% {
%   "maskSets": [
%     {
%       "maskFile": "mask1.tif",
%       "motherImageFile": "sample1.h5",
%       "maskOriginalLocation": "C:\\data\\mask1.tif",
%       "motherOriginalLocation": "C:\\data\\sample1.h5",
%       "type": "SAXS",
%       "timestamp": "2026-03-13 16:20:00"
%     }
%     {
%       "maskFile": "mask2.tif",
%       "motherImageFile": "sample1.h5",
%       "maskOriginalLocation": "C:\\data\\mask2.tif",
%       "motherOriginalLocation": "C:\\data\\sample1.h5",
%       "type": "SAXS",
%       "timestamp": "2026-03-13 16:20:00"
%     }
%   ],
%   "calibrationSets": [
%     {
%       "imageFile": "calib_img1.tif",
%       "setupFile": "calib_setup1.txt",
%       "imageOriginalLocation": "C:\\data\\calib_img1.tif",
%       "setupOriginalLocation": "C:\\data\\calib_setup1.txt",
%       "type": "SAXS",
%       "timestamp": "2026-03-13 16:20:05"
%     }
%   ]
% }


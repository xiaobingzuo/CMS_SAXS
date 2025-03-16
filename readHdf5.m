function h5data = readHdf5(h5file, beamline, choice, h5_keys)
% h5data = readHdf5(h5file, choice, h5_keys)
%   h5file: full name of hdf5 file
%   choice: 1: data file: image data and metaInfo data; 
%           2: master file: multiple data read, e.g., flatfield, mask, dataTime, etc
%   h5_keys: in Map structure: Key:value; i.e., {filed_name} & {path}
%     3/7/22   h5_keys is suppressed by "chioce". 


if nargin < 4
    if nargin < 3
        choice = 1;
    end
    h5_keys = eiger_preH5keys('data', beamline);
end

if choice == 1 
    h5_keys = eiger_preH5keys('data', beamline);
%     if ~isKey(h5_keys, 'data')
%         h5_keys = eiger_preH5keys();
%     end
    %h5data.data = h5read(h5file, h5_keys('data'));
    %h5data.phd  = h5read(h5file, h5_keys('It_phd'));
    keys = h5_keys.keys;
    for kk=1:h5_keys.Count
        key = keys{kk};
        cmd = sprintf('h5data.%s = h5read(h5file, h5_keys(''%s''));',key, key);
        try
            eval(cmd);
        catch
            warning('Field %s doesnot exist!', key);
        end
    end
end    
    
if choice == 2
    h5_keys = eiger_preH5keys('master', beamline);
%     if ~(isKey(h5_keys, 'flatfield') & isKey(h5_keys, 'mask'))
%          h5_keys = eiger_preH5keys();
%     end
    h5data.flatfield = h5read(h5file, h5_keys('flatfield'));
    h5data.mask = h5read(h5file, h5_keys('mask'));
    h5data.dataTime = h5read(h5file, h5_keys('dataTime'));
end


    function h5_keys=eiger_preH5keys(h5Data, beamline)
        h5_keys = containers.Map('KeyType','char','ValueType','char');
        
        switch h5Data
            case 'data'
                % data file        
                h5_keys('data') = '/entry/data/data';
                %%% metainfo keys
                %h5_keys('It_phd') = '/entry/Metadata/It_phd';
                %h5_keys('IC1_phd') = '/entry/Metadata/IC1_phd';
                
                %keys in MetaInfo
                if contains(beamline, '12-ID-B')
                    metaInfoKeys = {'AbsIntCoeff', 'AbsInt_Standard', 'Q_Standard', 'Beam_x_pixel', 'Beam_y_pixel', 'monoE','EnergyThres1','ExposureTime', 'Wavelength', 'Sample_Time', 'SDD', 'hexH', 'hexV', 'UserName', 'GUPNumber', 'ESAFNumber' };
                elseif contains(beamline, '12-ID-C') 
                    metaInfoKeys = {'It_phd', 'IC1_phd', 'AbsIntCoeff', 'AbsInt_Standard', 'Q_Standard', 'Beam_x_pixel', 'Beam_y_pixel', 'monoE', 'EnergyThres1','ExposureTime', 'Sample_Time', 'SDD', 'hexH', 'hexV', 'UserName', 'GUPNumber', 'ESAFNumber' };
                else
                    metaInfoKeys = {};
                end
                for kk =1:numel(metaInfoKeys)
                    key = metaInfoKeys{kk};
                    h5_keys(key) = sprintf('/entry/Metadata/%s',key);
                end
                
                %keys in MetaInfo
                if contains(beamline, '12-ID-B')
                    sampleInfoKeys = {'name', 'Io_flux', 'It_flux', 'sample_theta', 'sample_x_sth', 'sample_y_sav', 'sample_y_stv', 'sample_z_samx', 'thickness', 'temperature'};
                elseif contains(beamline, '12-ID-C')
                    sampleInfoKeys = {'name', 'Io_flux', 'It_flux', 'sample_theta', 'sample_x_sth', 'sample_y_sav', 'sample_y_stv',  'thickness', 'temperature'};
                else
                    sampleInfoKeys = {};
                end
                for kk =1:numel(sampleInfoKeys)
                    key = sampleInfoKeys{kk};
                    h5_keys(key) = sprintf('/entry/sample/%s',key);
                end
            case 'master'
                % master file
                detStr = '/entry/instrument/detector';
                detDescription = [detStr '/description'];
                detSDD = [detStr '/detector_distance'];
                h5_keys('detSDD') = detSDD;
                h5_keys('detDescription') = detDescription; 

                detSpecStr = '/entry/instrument/detector/detectorSpecific';
                flatfield = [detSpecStr '/flatfield'];
                dataTime = [detSpecStr '/data_collection_date'];
                h5_keys('flatfield') = '/entry/instrument/detector/detectorSpecific/flatfield';
                h5_keys('mask') = '/entry/instrument/detector/detectorSpecific/pixel_mask';
                h5_keys('dataTime') = dataTime;
            otherwise
                fprintf('%s: No such type in hdf5!\n', h5Data);
                
        end
    end


end


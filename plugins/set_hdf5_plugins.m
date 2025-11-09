function outputArg = set_hdf5_plugins()
% Set hdf5 plugins
% The hdf5 plugin dll files were obtained from The HDF Group.
% https://www.hdfgroup.org/solutions/hdf5/

inputArg = -1;

if isfolder(getenv("HDF5_PLUGIN_PATH"))
    inputArg = 200;
else
    fullPath = mfilename('fullpath');    
    [folderPath, ~, ~] = fileparts(fullPath);

    disp(['The folder name of the running command is: ', folderPath]);
    hdf5Str1 = sprintf('setx HDF5_PLUGIN_PATH "%s"', folderPath);
    inputArg = system(hdf5Str1);

    % hdf5Str2 = sprintf('setx HDF5_DIR "%s"', folderPath);
    % hdf5Str3 = sprintf('setx LD_LIBRARY_PATH "%s"', folderPath);    
    %system(hdf5Str2);
    %system(hdf5Str3);
end

if (inputArg ==200)
    fprintf("HDF5_PLUGIN_PATH already exists! Should ready to run.\n");
elseif (inputArg == 0)
    fprintf("HDF5 PLUGIN is set successfully!\n");
else
    fprintf("HDF5 PLUGIN is NOT set successfully!\n");
end

end
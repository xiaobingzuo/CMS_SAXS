function h5conv2tiff2(filename)
entry = '/entry/data/data';
R = h5read(filename, entry);
R = R';
[filepath,name,~] = fileparts(filename);
filename = fullfile(filepath, [name, '.tif']);
imwritetiff2(R, filename, 32)
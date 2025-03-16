function h5conv2tiff(filename)
entry = '/entry/data/data';
R = h5read(filename, entry);
R = R';
[filepath,name,~] = fileparts(filename);
filename = fullfile([name, '.tif']);
imwritetiff(R, filename, 32)
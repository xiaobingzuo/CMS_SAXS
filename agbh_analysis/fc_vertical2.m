function [p, m, g] = fc_vertical2(a, beamcenter, thres)
% [p, g] = fc_vertical2(img, beamcenter, thres)
% img : input diffraction image;
% beamcenter : beamcenter (x, y)
% thres : intensity threshold to define peaks.
% 
% p : 1D cell array, where each cell contains (x, y, d-spacing, th) positions
%   of diffraction peaks with the same d-spacing.
% m : 1D cell array, where each cell contains (x, y, d-spacing, th) positions
%   of diffraction peaks with the same azimuthal angle.
% g : 2D cell array, number of groups x number of angular sector
%   , of which each cell contains (x, y, d-spacing, th) positions
%   of diffraction peaks with the same d-spacing.
% Byeongdu Lee
% 3/26/2017

if nargin < 3
    thres = 1000;
end
x = 1:size(a, 1);
y = 1:size(a, 2);
[X,Y]=meshgrid(y,x);
t = a>thres;
x = X(t);
y = Y(t);
img = a(t);
%y = size(a,2)-y;
[th, d] = cart2pol(x-beamcenter(1), y-beamcenter(2));

% Sectoring
sec = linspace(-pi, pi, 180);
[nd, angbin] = histc(th, sec);
sec = find(nd > 0);


% grouping the peaks within similar d-spacings.
threscounts = 10;
g = {};
[n, BIN] = histc(d, 1:max(d));
nonzeros = n > threscounts;
indx = 2;
diff_nonzeros = diff(nonzeros);
indx_start = 0;
groupcount = 1;
newg = zeros(size(BIN));
while (indx < numel(nonzeros))
    if (diff_nonzeros(indx-1) == 1)
        indx_start = indx;
        newg = zeros(size(BIN));
    end
    if (indx_start > 0) && (diff_nonzeros(indx-1) == 0)
        newg = newg | BIN==indx;
    end
    if (indx_start > 0) && (diff_nonzeros(indx-1) == -1)
        %indx_end = indx-1;
        %g(groupcount) = indx_start:indx_end;
        for k=1:numel(sec)
            angsec = sec(k);
            g{groupcount, k} = find(newg & (angbin == angsec));
        end
        indx_start = 0; 
        groupcount = groupcount + 1;
    end
    indx = indx + 1;
end

% picking up only 1% intensity in each sector.
%figure(3)
for i = 1:size(g, 1)
%    c = rand(1, 3);
    p{i} = [];
    for j = 1:size(g, 2)
        indxs = g{i, j};
        if ~isempty(indxs)
            intv = img(indxs);
            [~, mvind] = max(intv);
            g{i,j} = [x(indxs(mvind(1))), y(indxs(mvind(1))),...
                d(indxs(mvind(1))), th(indxs(mvind(1)))];
            p{i} = [p{i}; g{i,j}];
            %mv = max(intv);
            %t = img(indxs) > mv*0.99;
            %g{i,j} = [x(indxs(t)), y(indxs(t))];
%            plot(g{i,j}(:,1), g{i,j}(:,2), 'o', 'color', c); hold on
        end
    end
end

for i = 1:size(g, 2)
    m{i} = [];
    for j = 1:size(g, 1)
        indxs = g{j, i};
        if ~isempty(indxs)
            m{i} = [m{i}; g{j,i}];
        else 
            m{i} = [m{i}; NaN, NaN, NaN, NaN];
        end
    end
end


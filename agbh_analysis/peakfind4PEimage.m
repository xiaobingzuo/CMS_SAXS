function g = peakfind4PEimage(a, center)
% [p, allpeaks] = peakfind4PEimage(a, center)
% a = imread('Pstandard_121.tif');
% [p, allpeaks] = peakfind4PEimage(flipud(a), [2172.7, 137.5])

siz = size(a);
if nargin < 2;
    center = [2172.7, 137.5];
end
%x = 1:550;
%% New ways;
%a = flipud(a);
% Search pixels that have counts higher than a threshold.
intensitythreshold = 2500;
t = a>intensitythreshold;
[I,J] = ind2sub(size(a),find(t==1));

% masking
t = I==3492;
I(t) = [];J(t) = [];

% histogram based on distances from beam center;
allpeaks = [I(:)-center(2), J(:)-center(1)];
p = sqrt(allpeaks(:,1).^2+allpeaks(:,2).^2);
[N,bin] = histc(p, 0:2:4000);
t = N > 30; % 

%grouping
g = [];k = 0;
m = diff(t);dt = [];p = {};
allpeaks = [];
% Grouping algorithm
% t is 1 when its count of pixels whose intensity is higher than the
% threshold in the range(which is 0:2:4000). The range containing less than
% 30 number of pixels does likely not have peaks but noise.
% So, we just need to group the ranges that have continously 1.
% To do this, t is differentiated and save to m.
% When m==0.5 and t == 0, new group starts and when m==-0.5 and t == 1 the
% group ends.
%
% After pixels are grouped, the intensities of the pixels in each group are
% investigated and the ones with higher intensities than the mean value of
% the group will remain selected.

for i=1:numel(m)
    if (t(i)==0) & (m(i)>0) 
        if ~isempty(dt)
            k = k+1;
            ind = zeros(size(bin));
            for i=1:numel(dt)
                ind = ind | bin == dt(i);
            end
            p{k} = [J(ind), I(ind)];
            
            % Selecting only pixels with higher intensity than the mean
            % value.
            ind2 = sub2ind(siz, I(ind), J(ind));
            allpeaks{k} = a(ind2);
            t2 = allpeaks{k} > mean(allpeaks{k})*1.0;
            p{k}(~t2, :) = [];
            allpeaks{k} = allpeaks{k}(t2);
            g(k).X = p{k}(:,1);
            g(k).Y = p{k}(:,2);
            g(k).val = allpeaks{k};
        end
        dt = [];
    end
    if (t(i)==1)
        dt = [dt, i];
    end
end
%p(1) = [];
%allpeaks(1) = [];
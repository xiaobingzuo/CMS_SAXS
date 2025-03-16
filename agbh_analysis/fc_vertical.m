function p = fc_vertical(a, step)
% p = fc_vertical(a, res, step, thres)
% find peaks of Agbh from WAXS on 12IDB...
% a is double(image);
% res = 5; % fpeak resolution
% step = 3;
% thres = 500;
%x = 1:size(a,1);
res = 2.5;

%x = 1:550;
x = 1:size(a, 1);
%% Peak searching
% For a vertical cut, background curve will be drawn with splinefitn.m
% Peaks then will be detected with fpeak.m for the background subtracted
% vertical cut.
fprintf('Peak searching started.\n');
peak = [];
for m=1:step:(size(a,2)-step+1)
    dt = sum(a(x, m:(m+step-1)), 2);
    pp1 = splinefitn(x, dt', 6, 2);back = ppval(pp1, x);
    dt = dt-back(:);
    %if sum(dt)>0
    t = fpeak(x, dt, res);
    %k = find((t(:,2) <= thres));
    %t(k,:) = [];
    k = (t(:,1) >= 0) & (t(:,1) <=1);
    t(k,:) = [];
    k = (t(:,1) >= 195) & (t(:,1) <=213);
    t(k,:) = [];
    k = (t(:,1) >= 407) & (t(:,1) <=425);
    t(k,:) = [];
    %plot(t(:,1), t(:,2), 'ro');hold on;
    peak = [peak;[ones(size(t(:,1)))*(m+step/2-1), t(:,1), t(:,2)]];
    %end
end
%N_count = m;
%% Remove peaks with intensity lower than the threshold.
% Threshold will be adjusted to keep the number of peaks less than 1800.
% If the number of peaks are larger than that, it will take very long time 
% to group them in the following step.
N_peak = numel(peak(:,1));
thres = 10;
while N_peak > 1800
    k = peak(:,3)>thres;
    peak = peak(k, :);
    N_peak = numel(peak(:,1));
    thres = thres + 10;
end
peak(:,3) = [];
fprintf('Peak searching done: %i peaks are found with threshhold intensity %i.\n', N_peak, thres);
assignin('base', 'allpeaks', peak)
fprintf('Sorting started.\n');
p = [];

%% Grouping peaks. 
% When two peaks is located at distance greater than the step and less than 
% 1.5 times of the step, they belong to the same group.

% Put first peak into the group #1.
p{1} = peak(1,:);peak(1,:) = [];
isPartofP = 1;
while ~isempty(peak) 
    while isPartofP % As long as one or more of peaks are found to be 
                    % a member(s) of existing groups. 
                    % Before next run, those assigned peaks will be removed.
        isPartofP = 0;
        removep = [];    
        for i=1:numel(peak(:,1))
            pd = peak(i,:);
            isPartofP = 0;
            for j=1:numel(p)
                t = p{j};
                d = sqrt((t(:,1)-pd(1)).^2+(t(:,2)-pd(2)).^2);
                tt = find((d >=step) & (d < 1.5*step), 1);
                if ~isempty(tt)
                    p{j} = [t;pd];
                    removep = [removep, i];
                    isPartofP = 1;
                end
            end
        end
        peak(removep,:) = []; % Peaks assigned to groups will be removed.
        if isempty(peak)
            break
        end
    end
    if isempty(peak)
        break
    end
    % Now if no particles are belong to existing groups, make a new group
    % with the first element of the peak.
    p{numel(p)+1} = peak(1,:);peak(1,:) = [];
    isPartofP = 1;
end

% Count the number of peaks in each group.
N = cellfun(@numel, p);
%k = N < N_count/5;
%p(k) = [];
%N = cellfun(@numel, p);
% Sort groups by the number of elements.
[~, bb] = sort(N, 'descend');
% Keep the top 10 groups
if numel(bb) > 10
    p = p(bb(1:10));
end
fprintf('Sorting Done.\n');
%N = zeros(size(p));
%for i=1:numel(p)
%    N(i) = numel(p{i})/2;
%end
assignin('base', 'p', p)
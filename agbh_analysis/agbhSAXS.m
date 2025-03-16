function [center, pixeldistance] = agbhSAXS(img, option)
% When a good guess is provided..
%[iniy, inix] = ginput(1);
%imgfilename = 'Sagbh_003_001.tif';
% option 1 or no option:
%   use auto analysis and if failed use clicking method.
% option 0
%   use only auto analysis
% option 2
%   use only the clicking method

LowerCntLimit = 250;
UpperCntLimit = 500;
img = double(img);
if nargin < 2
    option = 1;
end
Rfactor = 0;

if option < 2
    try
        [xc, yc, possibleR] = agbhImageprocess(img, LowerCntLimit, UpperCntLimit);
        Rfactor = 1;
    catch
        Rfactor = 0;
    end

    if Rfactor > 0
        img(img<0) = NaN;
        [inix, iniy, R, Rfactor] = radialcut_circlefit(img, xc, yc, possibleR);
    end
end

if (Rfactor < 5) && (option > 0);
    disp('Failed to find the center')
    %disp('Click the possible center')
    s = evalin('base', 'SAXSimageviewerhandle');
    %figure(s);
    %pnt = ginput(1);
    %xc = pnt(1); yc = pnt(2);
    disp('Pick 10 points on the first order peak at various directions')
    figure(s);
    pnt = ginput(10);
    xp = pnt(:,1);yp = pnt(:,2);
    [yc,xc,possibleR] = circfit(xp,yp);
    img(img<0) = NaN;
    [inix, iniy, R, ~] = radialcut_circlefit(img, xc, yc, possibleR);
end    
center = [iniy, inix];
pixeldistance = R;
fprintf('Refining xc = %0.2f, yc = %0.2f, 1st peak at %0.2f\n', inix, iniy, pixeldistance);


function [xc, yc, possibleR] = agbhImageprocess(a, LowerCntLimit, UpperCntLimit)
%a = imread(imgfilename);
%tic
[xdim, ~] = size(a);
if xdim==1679
    a = a(1:800, :);
    xdim = 800;
end
ress = [];
while isempty(ress)
    [ress, inix0, iniy0, t] = findcircle(a, LowerCntLimit, UpperCntLimit, xdim);
    LowerCntLimit = LowerCntLimit/3;
    UpperCntLimit = UpperCntLimit/3;
    if LowerCntLimit < 2
        break
    end
end

%findcircle
%if isempty(ress)
%    LowerCntLimit = LowerCntLimit/3;
%    UpperCntLimit = UpperCntLimit/3;
%    findcircle
%end

x = mod(t, xdim);
y = (t-x)/xdim + 1;
peak = [x, y];
Rt = sqrt((peak(:,1)-inix0).^2 + (peak(:,2)-iniy0).^2);

[n,x]=hist(Rt,fix(max(Rt)));
m = fpeak(x,n, 20);
m(m(:,2) < 50, :) = [];
%n(n<50) = 0;
%m = fpeak(x,n, 40);
Rpop = m(end, 1);
indx = find(x==Rpop);
indp = indx;
indm = indx;
peakindex = m(:,2)==max(m(:,2));
possibleR = m(peakindex, 1);
while 1
    indp = indp + 1;
    if indp > numel(n)
        indp = indp - 1;
        break
    end
    if n(indp) == 0
        break
    end
    if (indp-indx > 200)
        disp('Move the lower discriminator up 1')
        break
    end
end
while 1
    indm = indm - 1;
    if n(indm) == 0
        break
    end
    if (indx-indm > 200)
        disp('Move the lower discriminator up 2')
        break
    end
end
ind = (Rt>x(indm)) & (Rt<x(indp));
pn = peak(ind, :);
[xc,yc] = circfit(pn(:,1), pn(:,2));
fprintf('First Guess: xc = %0.2f, yc = %0.2f\n', (iniy0+1), (inix0+1));

% make a horizontal linecut and find peaks.
% assuming that the image has a ring and at least the first order peak of
% the ring shows up on the images at the two opposite sides of the direct beam.
resol = 20;xax = 1:size(a, 2);numpeak = 0;
ROI = fix(xc)-2:fix(xc)+2; % 5 rows will be added to improve statistics.
while ((numpeak < 2) || (numpeak > 12))
    t = fpeak(xax, sum(a(ROI, :)), resol);
    if ~isempty(t)
        numpeak = numel(t(:,1));
    else
        numpeak = 0;
    end
    if numpeak < 2
        resol = resol - 10;
    else
        resol = resol + 10;
    end
    if (resol < 10) || (resol > 100)
        break
    end
end
[y, i] = sort(t(:,end), 1, 'descend');
x = t(i,1);
% Let's validate the first two peaks are the two peaks we are looking for.
% Intensities of the first two peaks should not be so much different.
% And, their distances from the center should not be different either. 

% First remove peaks so close to the beam position, which may be due to the
% beamstop..
kk = find(abs(x-yc) < 20);
x(kk) = []; y(kk) = [];

isok = 0;
while ~isok
    if (y(1) < 3*y(2)) && (abs(sum(x(1:2))-2*yc) < 30)
        possibleR = sum(abs(yc-x(1:2)))/2;
%        isok = 1;
        break
    else
        isok = 0;
    end
    if (numel(x) < 3)
        break
    end
    y(1) = [];
    x(1) = [];
end

%if abs(sum(x(1:2))-2*yc) < 10
    % provide new possibleR.
%    possibleR = sum(abs(yc-x(1:2)))/2;
%end

fprintf('First Fit: xc = %0.2f, yc = %0.2f, R=%0.2f\n', (yc), (xc), possibleR);
%toc

function [ress, inix, iniy, tx] = findcircle(a, LowerCntLimit, UpperCntLimit, xdim)
    b=ones(size(a));
    b(a<LowerCntLimit) = 0; b(a>=UpperCntLimit) = 0;
    tx = find(b==1);

    test = tx;
    while 1
        if numel(test)<1
            disp('change high and low limit')
            inix = 1;
            iniy = 1;
            break
        end
        pnn = maskfindpeak(b, test(1),10);
        if numel(pnn) > 2000
            x0 = mod(pnn, xdim);
            y0 = (pnn-x0)/xdim + 1;
            p = [x0, y0];
            [inix,iniy] = circfit(p(:,1), p(:,2));
            break
        end
        test = setdiff(test, pnn);
    end
    ress=test;
end



end
end


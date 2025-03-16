function [xc,yc,indx] = findPeakCenter(p1, p2, img)
%[xc,yc] = findPeak(p1, p2, img)
% p1=[x1,y1]; p2 = [x2,y2]; all whole numbers
%width: band width in pixel#
%img: 2D image data
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

figH = evalin('base', 'SAXSimageviewerhandle');handles = guidata(figH);
saxs = getgihandle;

indx=1;
xc = 1;
yc=1;

%img=flipud(img);
figure; imagesc(img);
[col, row]=size(img);
p1x=p1(1);
p1y=p1(2);
p2x=p2(1);
p2y=p2(2);
width = round(sqrt(sum((p1-p2).^2)));
if max(p1x, p2x)>col | min(p1x, p2x)<1 | max(p1y, p2y)>row | min(p1y, p2y)<1
    indx = 0;
    return;
end

x = round(linspace(p1x,p2x,width));
y = round(linspace(p1y,p2y,width));
z=[];
for kk=1:width
    %z=[z img(x(kk),y(kk))];
    z=[z img(y(kk),x(kk))];
end
zx=1:width;
tmph1 = line(x,y, 'parent', handles.ImageAxes, 'color', 'r');

figure(501); plot(zx,z,'ob');
%f = ezfit(zx, z, 'gauss');
f = ezfit(zx, z, 'y(x) = a*exp(-(x-x0)^2/(2*s^2))+c0');
%       param: {'a'  'c0'  's'  'x0'}
%           m: [1.2522e+03 1.0595e+02 2.2074e+00 2.8834e+01]
showfit(f, 'fitcolor', 'red', 'fitlinewidth', 2); 
if f.r < 0.90
    indx = 0;
else
    zc = round(f.m(end));
    if zc*(width-zc)>=0
        xc = x(zc); yc = y(zc);
        indx = 1;
    else
        indx = 0;
    end
end
end


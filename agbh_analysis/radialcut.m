function [R, Iq, x, y] = radialcut(img, pnt, tiltang, bandwidth)
% [R, Iq] = radialcut(img, pnt, tiltang)
if nargin < 4
    bandwidth = 1;
end
[mx, my] = size(img);
diag=round(sqrt(mx^2+my^2));
R = linspace(0, diag,(diag+1)*5);
phi = deg2rad(tiltang)*ones(size(R));
[x,y] = pol2cart(phi,R);
x = x+pnt(1);
y = y+pnt(2);
t1 = find((x>my)|(x<0));
t2 = find((y>mx)|(y<0));
t = [t1(:);t2(:)];
t = unique(t);
x(t) = [];x = x(:);
y(t) = [];y = y(:);
R(t) = [];R = R(:);

if bandwidth>1
    band = -bandwidth:bandwidth;
    xt = cos(tiltang*pi/180+pi/2);
    yt = sin(tiltang*pi/180+pi/2);
    x2 = repmat(x, 1, numel(band)) + repmat(band*xt, numel(x), 1);
    y2 = repmat(y, 1, numel(band)) + repmat(band*yt, numel(x), 1);
    %x2 = repmat(x, 1, numel(band)) + repmat(band*cos(tiltang*pi/180+pi), numel(x), 1);
    %y2 = repmat(y, 1, numel(band)) + repmat(band*sin(tiltang*pi/180+pi), numel(x), 1);
else
    x2 = x;
    y2 = y;
end
Iq = sum(interp2(img, x2, y2),2);
Iq(isnan(Iq)) = 0;

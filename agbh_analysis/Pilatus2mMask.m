function res = Pilatus2mMask(pntX, pntY)
% function res = Pilatus2mMask(pntX, pntY)
% pnt is a point by coordinate.
% If pnt is on a nonactive area, return 1 else 0.
% Example:
%   Pilatus2mMask([], 50)
%
if isempty(pntX)
    pntX = 1;
end
if isempty(pntY)
    pntY = 1;
end

Pilatus2m.Xactive = 194;
Pilatus2m.Xnonactive = 18;
Pilatus2m.Xnum = 7;
Pilatus2m.Xperiod = Pilatus2m.Xactive+Pilatus2m.Xnonactive;

Pilatus2m.Yactive = 486;
Pilatus2m.Ynonactive = 8;
Pilatus2m.Ynum = 2;
Pilatus2m.Yperiod = Pilatus2m.Yactive+Pilatus2m.Ynonactive;

Xbad = [];
for i=1:Pilatus2m.Xnum
    x0 = Pilatus2m.Xperiod*(i-1) + Pilatus2m.Xactive+1;
    y0 = x0 + Pilatus2m.Xnonactive;
    Xbad = [Xbad; [x0, y0]];
end
Ybad = [];
for i=1:Pilatus2m.Ynum
    x0 = Pilatus2m.Yperiod*(i-1) + Pilatus2m.Yactive+1;
    y0 = x0 + Pilatus2m.Ynonactive;
    Ybad = [Ybad; [x0, y0]];
end
Xbad = [Xbad(:,1)-1, Xbad(:,2)+1];
Ybad = [Ybad(:,1)-1, Ybad(:,2)+1];

InorNot = [(Xbad(:,1) <= pntY) & (Xbad(:,2) >= pntY);
    (Ybad(:,1) <= pntX) & (Ybad(:,2) >= pntX)];
if sum(InorNot) > 0
    res = 1;
else
    res = 0;
end


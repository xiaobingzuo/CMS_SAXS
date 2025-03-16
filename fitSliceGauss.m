function [found, peak] = fitSliceGauss(pt1, pt2, bandWidth, img, imgSize, fitR, mask)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 7
    mask = [];
end
found = 0;
peak = [];
%fitR = 0.9;

if bandWidth<4
    bandWidth =20;
end
xL =round(pt1(1)); yL =round(pt1(2));
xH =round(pt2(1)); yH =round(pt2(2));
 
if logical(size(mask) == imgSize) |  logical(isempty(mask))  
    if min(xL,xH)>0 && max(xL, xH)<imgSize(2) && min(yL,yH)>0 && max(yL, yH)< imgSize(1)
        x=round(linspace(xL,xH, bandWidth));
        y=round(linspace(yL,yH, bandWidth));
        %[xL xH yL yH]
        xySlice = [x;y]';
        if ~checkInMask(xySlice, mask) % this slice not in mask, continue to find peak
            zx=1:bandWidth;
            zy=[];
            for nn=1:bandWidth
                zy=[zy img(y(nn),x(nn)) ];
            end   
            
            ezfitStr = sprintf("f=ezfit(zx,zy,'zy(zx)=a*exp(-(zx-x0)^2/(2*s^2))+b; x0 =%f; a=%f; s=%f');", 0.5*bandWidth, max(zy), 5);
            eval(ezfitStr);
            %figure(555);clf; plot(zx,zy,'ro'); showfit(f);pause(0.5);
            %f.r
            xCnt = f.m(4);
            if f.r>fitR && (xCnt*(xCnt-bandWidth)<=0)
                found = 1;
                %peak = [x(round(xCnt));y(round(xCnt))];
                px=polyfit(zx,x,1); py=polyfit(zx,y,1);
                xt = polyval(px, xCnt); yt = polyval(py, xCnt);
                peak = [xt; yt];
            end
        end
        
    end
end
end


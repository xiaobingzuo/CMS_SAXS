function [ring, ring1, ring2] = calculateQringsYaw(qValue, pixelSize, beamXY, SDD, xeng, nPts, bandWidth, yaw, det)
% [ring, ring1, ring2] = calculateQrings(qValue, pixelSize, beamXY, SDD, xeng, nPts, bandWidth)
%
%   Detailed explanation goes here
if nargin ==8
    det = [487 619]; % pilatus300k yrow=619, xcol=487
end
yRow=det(2);
xCol=det(1);
xc=beamXY(1);
yc=beamXY(2);
SDDrel = SDD / pixelSize;
theta2 = q2rad2(qValue,xeng);
yawRad = yaw / 180*pi;

ring  = [];
ring1 = [];
ring2 = [];
% xe =[]; ye =[];
% xe1=[]; ye1=[];
% xe2=[]; ye2=[];
extraPts = 10;
%x = unique(round(linspace(1-extraPts,xCol+extraPts, nPts)))';
x = unique((linspace(1-extraPts,xCol+extraPts, nPts)))';

 for kk=1:length(x)  % lower half rings
     xval = x(kk);
     t = tan(theta2);
     sy = sin(yawRad);
     cy = cos(yawRad);
     
     A = (t *sy)^2 - (cy)^2;
     B = 2. * SDDrel * sy * t*t;
     C = (t*SDDrel)^2 -(xval-xc)^2;
     yRoots = roots([A B C]);
     if isreal(yRoots)
         yRs = yRoots(yRoots>0 & yRoots<yRow); % 1 solution
         
         % middle ring
         yval = yc + min(yRoots);
         ring = [ring; [xval yval]];

        % smaller ring
        phiRad = atan((yval -yc)/(xval-xc));
        dX = 0.5* bandWidth * (cos(phiRad));
        dY = 0.5* bandWidth * (sin(phiRad));    
        ring1 = [ring1; [xval+dX yval+dY]];
        
        % larger ring
        ring2 = [ring2;  [xval-dX yval-dY]];

     end
 end
 
 for kk=length(x):-1:1  % upper half rings
     xval = x(kk);
     t = tan(theta2);
     sy = sin(yawRad);
     cy = cos(yawRad);
     
     A = (t *sy)^2 - (cy)^2;
     B = 2. * SDDrel * sy * t*t;
     C = (t*SDDrel)^2 -(xval-xc)^2;
     yRoots = roots([A B C]);
     if isreal(yRoots)
         yRs = yRoots(yRoots>0 & yRoots<yRow); % 1 solution
         
         % middle ring
         yval = yc + max(yRoots);
         ring = [ring; [xval yval]];

        % smaller ring
        phiRad = atan((yval -yc)/(xval-xc));
        dX = 0.5* bandWidth * (cos(phiRad));
        dY = 0.5* bandWidth * (sin(phiRad));    
        ring1 = [ring1; [xval-dX yval-dY]];
        
        % larger ring
        ring2 = [ring2;  [xval+dX yval+dY]];

     end
 end

% y = unique(round(linspace(1,yRow, nPts)))';
% for kk=1:length(y)
%     yval = y(kk);
%     delta = (tan(theta2) * (SDDrel - (yval -yc)*sin(yawRad)))^2 -((yval -yc)*cos(yawRad))^2;
%     if delta >=0
%         % middle ring
%         xval1 = xc - sqrt(delta);
%         xval2 = xc + sqrt(delta);
%         ring = [ring; [xval1 yval]; [xval2 yval]];
%         
%         % smaller ring
%         phiRad = atan((yval -yc)/(xval1-xc));
%         dX = 0.5* bandWidth * sin(phiRad);
%         dY = 0.5* bandWidth * cos(phiRad);    
%         ring1 = [ring1; [xval1-dX yval-dY]; [xval1+dX yval+dY]];
%         
%         % larger ring
%         phiRad = atan((yval -yc)/(xval2-xc));
%         dX = 0.5* bandWidth * sin(phiRad);
%         dY = 0.5* bandWidth * cos(phiRad);    
%         ring2 = [ring2; [xval2-dX yval-dY]; [xval2+dX yval+dY]];
%                 
%     end
% 
% end
%  for kk=1:length(x)
%      xval = x(kk);
%      t = tan(theta2);
%      sy = sin(yawRad);
%      cy = cos(yawRad);
%      
%      A = (t *sy)^2 - (cy)^2;
%      B = -2. * SDDrel * sy * t*t;
%      C = (t*SDDrel)^2 -(xval-xc)^2;
%      yRoots = roots([A B C]);
%      if isreal(yRoots)
%          yRs = yRoots(yRoots>0 & yRoots<yRow); % 1 solution
%          
%          % middle ring
%          yval1 = yc + yRoots(1);
%          yval2 = yc + yRoots(2);
%          ring = [ring; [xval yval1]; [xval yval2]];
% 
%         % smaller ring
%         phiRad = atan((yval1 -yc)/(xval-xc));
%         dX = 0.5* bandWidth * sin(phiRad);
%         dY = 0.5* bandWidth * cos(phiRad);    
%         ring1 = [ring1; [xval-dX yval1-dY]; [xval+dX yval1+dY]];
%         
%         % larger ring
%         phiRad = atan((yval2 -yc)/(xval-xc));
%         dX = 0.5* bandWidth * sin(phiRad);
%         dY = 0.5* bandWidth * cos(phiRad);    
%         ring2 = [ring2; [xval-dX yval2-dY]; [xval+dX yval2+dY]];
%      end
% end
% ring  = sortrows(ring, 2);
% ring1 = sortrows(ring1, 2);
% ring2 = sortrows(ring2, 2);

% th = linspace(0, 2*pi, nPts);
% th = th';
% Re = round(tan(q2rad2(qValue,xeng)) * (SDD / pixelSize)); 
% xe = Re*cos(th)+ beamXY(1); 
% ye = Re*sin(th)+ beamXY(2);
% 
% ring = [xe ye];
% 
% Re1 = Re - bandWidth / 2;
% Re2 = Re + bandWidth / 2;
% 
% xe1 = Re1*cos(th)+ beamXY(1); 
% ye1 = Re1*sin(th)+ beamXY(2);
% 
% xe2 = Re2*cos(th)+ beamXY(1); 
% ye2 = Re2*sin(th)+ beamXY(2);
% 
% ring1 = [xe1 ye1];
% ring2 = [xe2 ye2];
end


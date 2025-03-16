function [ring, ring1, ring2] = calculateQrings(qValue, pixelSize, beamXY, SDD, xeng, nPts, bandWidth, yaw)
% [ring, ring1, ring2] = calculateQrings(qValue, pixelSize, beamXY, SDD, xeng, nPts, bandWidth)
%
%   Detailed explanation goes here
if nargin ==7
    yaw = 0;
end
th = linspace(0, 2*pi, nPts);
th = th';
Re = round(tan(q2rad2(qValue,xeng)) * (SDD / pixelSize)); 
xe = Re*cos(th)+ beamXY(1); 
ye = Re*sin(th)+ beamXY(2);

ring = [xe ye];

Re1 = Re - bandWidth / 2;
Re2 = Re + bandWidth / 2;

xe1 = Re1*cos(th)+ beamXY(1); 
ye1 = Re1*sin(th)+ beamXY(2);

xe2 = Re2*cos(th)+ beamXY(1); 
ye2 = Re2*sin(th)+ beamXY(2);

ring1 = [xe1 ye1];
ring2 = [xe2 ye2];
end


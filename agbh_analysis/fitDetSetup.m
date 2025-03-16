function cv = fitDetSetup(p, XY, qval, otherFixedParameters)
% Fitting function for Detector from peak positions found.
% cv = cv = fitDetSetup(p, XY, qval);
% p : Setup Parameters
%    centerX
%    centerY
%    sdd
%    Pitch
%    Yaw
%    Roll
% XY : pixel coordinates, a 1xm cell, where m is the total number of q
% values.
% qval : q values for XYs. 1xm array
% otherFixedParameters
%   otherFixedParameters.pixelsize
%   otherFixedParameters.waveln


% Input of pixel2q:
%     pix = varargin{1};
%     center = varargin{2};
%     sdd = varargin{3};
%     pixelsize = varargin{4};
%     detector_tiltangle = varargin{5};
%     lambda = varargin{6};
% 
center = [p(1), p(2)];
sdd = p(3);
detangle = [p(4), p(5), p(6)];

psize = otherFixedParameters.pixelsize;
wl = otherFixedParameters.waveln;
err = 0;
dt = [];
qv = [];
for i=1:numel(XY)
    pix = XY{i};
    %q = pixel2q(pix, center, sdd, psize, detangle, wl);
    %pix = bsxfun(@minus,pix,center);
    if isempty(pix)
        continue
    end
    if size(pix, 2) == 4
        pix = pix(:, 1:2);
    end
    pix = pix - repmat(center, length(pix(:,1)), 1);
    qv0 = pixel2qv0(pix, wl, sdd, psize, detangle);
    q = sqrt(sum(qv0.*qv0, 2));
    %err = err + sum((q-qval(i)).^2);
    dt = [dt;q];
    qv = [qv;qval(i)*ones(size(q))];
end
cv = chi_squared(qv, dt, 5);

%cv = chi_squared(yarray, Iqarray, 5);

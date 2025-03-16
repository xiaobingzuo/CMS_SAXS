function [inix, iniy, R, Rfactor] = radialcut_circlefit(a, xc, yc, pixeldistance)
% number of linecuts.
Numcut = 20;Hdir = 4;
bandwidth = 1; % linecut bandwidth

% horizontal.
%BeamstopRadius = 30;
%ang = -20:20:210;
%Xpos = [];    
% pixeldistance is required
ppos = [];
ypos = [];
figure;
ang = linspace(-10, 190, Numcut);
fprintf('fitting curves.')
Rfactor = 0;
for i=1:numel(ang)
    [R, Iq, xp, yp] = radialcut(a, [yc,xc], ang(i), bandwidth);
    indexout = R>pixeldistance*1.5;
    R(indexout) = [];
    Iq(indexout) = [];
    xp(indexout) = [];
    yp(indexout) = [];
    indexout = R<pixeldistance*0.75;
    R(indexout) = [];
    Iq(indexout) = [];
    xp(indexout) = [];
    yp(indexout) = [];
    
    %Iq = Iq.*R(:).^2;
    maxIq = max(Iq);
    k = find(Iq < maxIq*0.10);
    Iq(k) = []; R(k) = [];xp(k) = [];yp(k) = [];
    Imax = max(Iq);
    fprintf('.')
    if isempty(Iq)
        continue;
    end
    if isempty(R)
        continue;
    end
    minv = min(Iq);
    if numel(minv) > 2;
        minv = minv(1);
    end
    if isnan(minv);
        minv = 0;
    end
    minv = 0;
    f = ezfit(R, Iq, sprintf('y = a*exp(-1/2*(x-c)^2/b^2)+d*x+e; a=%f; b=2; c=%f; d=-2; e=%f', Imax/2, pixeldistance, minv));
    subplot(Hdir,Numcut/Hdir,i);
    if f.m(3) < 0;
        f.r = 0;
    end
    Rfactor = Rfactor + f.r;
    if f.r>0.9
        plot(R, Iq, 'bo');
    else
        plot(R, Iq, 'ro');
    end
    showfit(f, 'dispfitlegend','off', 'dispeqmode', 'off', 'dispeqboxmode', 'off')
    %showfit(f, 'dispfitlegend','off', 'dispeqmode', 'off')
    legend off;
    if f.r > 0.9
        xpoint = interp1(R, xp, f.m(3));
        ypoint = interp1(R, yp, f.m(3));
        %ppos = [ppos, f.m(3)];
        %ypos = [ypos, f.m(1)];
        ppos = [ppos, xpoint];
        ypos = [ypos, ypoint];
    end
    %Xpos{i} = ppos;
    %Ypos{i} = ypos;
end
fprintf('done\n')
k = isnan(ppos);
ppos(k) = [];
ypos(k) = [];
k = isnan(ypos);
ppos(k) = [];
ypos(k) = [];
[inix,iniy, R] = circfit(ppos, ypos);

%center = [inix, iniy];
pixeldistance = R;
fprintf('Refining xc = %0.2f, yc = %0.2f, 1st peak at %0.2f\n', yc, xc, pixeldistance);

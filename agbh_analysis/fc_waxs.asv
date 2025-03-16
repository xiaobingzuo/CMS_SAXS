a = imread('Wagbh_003_001.tif');
a = double(a);
[sx,sy] = size(a);
k=5;
Hcut = a(k, :);
m = fpeak(x, Hcut, 80); 
meanback = mean(Hcut);
m(find(m(:,2) < meanback*1.5), :) = [];
dist2center = (m(:,1) - sy/2);
k = dist2center > 0;
pm = m(k,:);
dt1 = pm(find(pm(:,1) == max(pm(:,1))),:);
k = dist2center <=0;
pn = m(k,:);
dt2 = pn(find(pn(:,1) == max(pn(:,1))),:);
dt = [dt1; dt2];

plot(x, Hcut, 'b', m(:,1), m(:,2), 'ro');
ppos = [];
for j=1:numel(dt(:,1))
    f = ezfit(x, Hcut, sprintf('y = a*exp(-1/2*(x-c)^2/b^2)+back; a=%f; b=2; c=%f; back = %f', dt(j,2), dt(j,1), meanback));
    ppos = [ppos, f.m(4)];
end

xc = mean(ppos);

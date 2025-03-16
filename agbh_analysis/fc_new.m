% When no ini point is provided..
tic
[xdim, ydim] = size(a);
b=ones(size(a));b(a<250) = 0; b(a>=500) = 0;
t = find(b==1);
i = 1;
xc = []; yc = [];
while numel(xc)<5
    pn = maskfindpeak(b, t(1),10);
    if numel(pn) > 2000
        x = mod(pn, xdim);
        y = (pn-x)/xdim + 1;
        p = [x, y];
        [xc(i),yc(i)] = circfit(p(:,1), p(:,2));
        i = i + 1;
    end
    t = setdiff(t, pn);
end
xc(outlier2(xc, 1, 1)) = [];
yc(outlier2(yc, 1, 1)) = [];

fprintf('First Guess: xc = %0.2f, yc = %0.2f\n', mean(xc), mean(yc));


%x = mod(t, xdim);
%y = (t-x)/xdim + 1;
%peak = [x, y];





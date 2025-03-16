function pn = maskfindpeak(maskimg, startingpnt, howfararound)

%[y, x] = ginput(1);
b = maskimg;
xdim = size(b,1);
ydim = size(b,2);
numstart = 2;
if numel(startingpnt) == 1
    x = mod(startingpnt, xdim);
    y = (startingpnt-x)/xdim + 1;
    startingpnt = [x, y];
    numstart = 1;
end
x = startingpnt(1);
y = startingpnt(2);
if nargin < 3
    howfararound = 2;
end

x = fix(x);y = fix(y);
if b(fix(x),fix(y))< 1
    disp('select brown region')
    return
end
%x0 = x;
%y0 = y;
%p = [x, y];
mask1 = b;
mask1(x,y) = 0.5;
%mask2 = mask1;
%psel =[];

shiftp = zeros((2*howfararound+1)^2-1, 2);
k = 1;
for i=-howfararound:1:howfararound
    for j=-howfararound:1:howfararound
        shiftp(k, 1) = i;
        shiftp(k, 2) = j;
        if (i~=0) || (j~=0)
            k = k+1;
        end
    end
end
m = 1;
nshiftp = numel(shiftp(:,1));

while (1)
    pn = [];
    %sprintf('Size of x data is %i at time %i', numel(x), m)
    m = m+1;
    xn = repmat(x, [nshiftp, 1]);
    yn = repmat(y, [nshiftp, 1]);
    spx = repmat(shiftp(:,1), [numel(x), 1]);
    spy = repmat(shiftp(:,2), [numel(x), 1]);
    pn = [pn;xn+spx, yn+spy];
    p = xdim*(pn(:,2)-1) + pn(:,1);
    p = unique(p);
    x = mod(p, xdim);
    y = (p-x)/xdim + 1;
    pn = [x, y];
    
    %for i=1:numel(x)
    %    p = [x(i)+shiftp(:,1), y(i)+shiftp(:,2)];
    %    pn = [pn;p];
    %end
    k = find((pn(:,1) < 1) | (pn(:,2) < 1));
    if ~isempty(k)
        pn(k,:) = [];
    end
    k = find((pn(:,1) > xdim) | (pn(:,2) > ydim));
    if ~isempty(k)
        pn(k,:) = [];
        %p(k) = [];
    end
    %valmask1 = [];
    p = xdim*(pn(:,2)-1) + pn(:,1);
    %for i=1:size(pn,1)
    %    valmask1(i) = mask1(pn(i,1), pn(i,2));
    %end
    valmask1 = mask1(p);
    k = find(valmask1 == 0.5);
    if ~isempty(k)
        %pn(k,:) = [];
        p(k) = [];
    end
    %valmask2 = [];
    %p = xdim*(pn(:,2)-1) + pn(:,1);
    %for i=1:size(pn,1)
    %    valmask2(i) = b(pn(i,1), pn(i,2));
    %end
    valmask2 = b(p);
    k = find(valmask2 < 1); %if pn is not on the masked area.
    if ~isempty(k)
        %pn(k,:) = [];
        p(k) = [];
    end
    %p = xdim*(pn(:,2)-1) + pn(:,1);
    
    mask1(p) = 0.5;
%    psel = [];
    p = unique(p);
    x = mod(p, xdim);
    y = (p-x)/xdim + 1;
    pn = [x, y];
    %if ~isempty(pn)
        %for i=1:numel(pn)/2
        %    mask1(pn(i,1), pn(i,2)) = 0.5;
        %end
    %    psel = [psel; pn];
    %end
    
    if isempty(pn)
        break
    end
    %x = psel(:,1);
    %y = psel(:,2);
    %psel = [];
        
    %mask2 = b;
    %mask2(p) = 0.5;
    % set(findobj(gca, 'type', 'image'), 'Cdata', mask1);
    %drawnow
end

pn = find(mask1 == 0.5);
if numstart == 2
    x = mod(pn, xdim);
    y = (pn-x)/xdim + 1;
    pn = [x, y];
end


        
    
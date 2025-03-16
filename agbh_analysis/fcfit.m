function ret = fcfit(x, peak)
ret = sum((peak(:,1)-x(1)).^2 + (peak(:,2)-x(2)).^2);

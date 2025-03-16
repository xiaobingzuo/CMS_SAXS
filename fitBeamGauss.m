function fitxy = fitBeamGauss( data, params, smth )
% fitxy = fitBeamGauss( data, params, smooth )
%   data = [x y] to fit with gauss+y0
%   params = [y0 a x0 sigma]
%   fit_fun: y(x) = y0 + a*exp(-(x-x0)^2/(2*sigma^2));
%   fitxy = [x y_fit]; y_fit = y0 + a*exp(-(x-x0)^2/(2*sigma^2));
%   smooth=1: need to smooth data first; 
%     smooth=0 or only two arguments: no smooth

if nargin==2
    smth =0 ;
else
    smth =1;
end

x = data(:,1);
y = data(:,2);
if smth ==1
    y = smooth(y);
end
y0 = params(1);
a  = params(2);
x0 = params(3);
sigma = params(4);
fun=sprintf('y(x) = y0 + a*exp(-(x-x0)^2/(2*sigma^2)); y0=%f; a=%f; x0=%f; sigma=%f', y0, a, x0, sigma);
f = ezfit(x, y, fun);

% fitting values
eval(sprintf('%s=%f',f.param{1},f.m(1)));
eval(sprintf('%s=%f',f.param{2},f.m(2)));
eval(sprintf('%s=%f',f.param{3},f.m(3)));
eval(sprintf('%s=%f',f.param{4},f.m(4)));

y_fit = y0 + a*exp(-(x-x0).^2/(2*sigma^2));
fitxy=[x y];
figure(998); plot(x,y, 'or', x,y_fit, '-b'); title(sprintf('beam FWHM: %.3f', 2.355*sigma));

%
% display fitting
figure(999); plot(x,y, 'or'); showfit(f)
sprintf('beam size FWHM 2.355*sigma : %.3f', 2.355*sigma)


end
function theta2 = q2rad2( q, eng )
% theta2 = q2angle2( q, eng )
%   q = 4*pi/wl *sin(theta)
%   return 2*theta

theta2=2.*asin(q./(4*pi/eng2wlen(eng)));
end
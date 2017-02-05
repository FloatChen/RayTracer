function [yout, xout]=rk4(y,dydx,x,h)
% [yout]=rk4_m(y,h,ellipseCon,fr)
% 4th-order Runge Kutta integration - see Numerical Recipies, p. 706.
% OR Wikipedia - Runga Kutta 4th Order
% y is vector [alpha lat lon]
% yout is vector [alpha lat lon] 

xout = x+h;

h2=h*0.5;
h6=h/6.0;

k1 = dydx;

k2 = getderivs(y+k1*h2);

k3 = getderivs(y+k2*h2);

k4 = getderivs(y+k3*h);

yout = y + h6*(k1 + 2*k2 + 2*k3 + k4);



end
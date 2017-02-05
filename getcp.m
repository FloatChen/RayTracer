function cp = getcp(f,h) 
% [cp, cg] = getc(f,h) 
%g = 9.81;
k = getk(f,h);
%cg = 0.5*(g*tanh(k*h)+g*(k*h)*(sech(k*h)^2))/sqrt(g*k*tanh(k*h));
%cp = sqrt(g*k*tanh(k*h))/k;
cp = 2*pi*f/k;

end
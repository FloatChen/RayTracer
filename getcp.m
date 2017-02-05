function cp = getcp(f,h) 
% [cp, cg] = getc(f,h) 

k = getk(f,-h);
cp = 2*pi*f/k;

end
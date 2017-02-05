function [ Yout, xout, hnext, maxcount ] = rkqc(Y,dydx,x,htry,eps,Yscale,hmin)
%[ Yout, xout, hnext] = rkqc(Y,dydx,x,htry,eps,yscal)


%Constants
safety = .9;
pshrink = -.25;
pgrow = -.20;
fcor = 1/15;
errcon = 6e-4;

%intialize
hnext = 0;
maxcount = 0;
while hnext == 0 && maxcount < 5
    maxcount = maxcount+1;
    h = htry;
    hh = .5*htry;
    
    % Short Step
    [Yshort, xshort] = rk4(Y,dydx,x,hh);
    dydx_2 = getderivs(Yshort);
    [Yshort, ~] = rk4(Yshort,dydx_2,xshort,hh);
    
    % Larg Step
    [Ylong, xout] = rk4(Y,dydx,x,h);
    
    % Evaluate accuracy
    Yerr = Ylong - Yshort;
    errmax = max(abs((Yshort-Ylong)./Yscale));
    
    % normalize error
    errmax = errmax/eps;
    
    if errmax > 1 %If too big make smaller and try again
        htry = safety*h*(errmax^pshrink);
        if htry <= hmin
            hnext = hmin;
        end            
    else    % Otherwise make step bigger and continue
        if errmax > errcon
            hnext = safety*h*errmax^pgrow;
        else
            hnext = 4*h; % Don't grow too fast
        end
    end
    
end

% Adjust to 5th order RungaKutta
%Yout = Yshort;
Yout = Yshort + Yerr*fcor;


end


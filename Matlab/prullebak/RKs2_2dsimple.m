function y = RKs2_2dsimple(xspan,h,alpha,y0,f)
global a
xn = xspan(1);
y = y0.';
sol = [0.1,0.1,0.1,0.1,0.1,0.1];
z = a*1i*h;
th = NRtheta(z^2,alpha,sqrt(3)/6);
while xn < xspan(2)-h/2
    sol = fsolve(@(Y) RK2system(Y,z,h,f,y(:,end),th),sol,optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16));
    xn = xn+h;
    y = [y,[sol(3);sol(6)]];
end
end
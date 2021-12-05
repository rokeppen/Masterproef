function y = RKs2_kepler(xspan,h,alpha,y0,f)
xn = xspan(1);
y = y0.';
sol = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
while xn < xspan(2)-h/2
    z=1i*h*(y(1,end)^2+y(2,end)^2)^(-3/2); %2 naar 4?
    th=NRtheta(z^2,alpha,sqrt(3)/6);
    sol = fsolve(@(Y) RK2system(Y,abs(z),h,f,y(:,end),th),sol,optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16));
    xn = xn+h;
    y = [y,[sol(3);sol(6);sol(9);sol(12)]];
end
end
function y = RKs2_simple(xspan,h,alpha,y0,f)
xn = xspan(1);
y = y0.';
sol = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
z = 1i*h;
th = NRtheta(z,alpha,sqrt(3)/6)
while xn < xspan(2)-h/2
    xn = xn+h;
    sol = fsolve(@(Y) RK2system(Y,z,h,f,y(:,end),th,xn),sol,optimoptions('fsolve','Display','Off'));
    y = [y,[sol(3);sol(6);sol(9);sol(12)]];
end
end
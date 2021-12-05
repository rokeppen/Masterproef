function y = NRtheta3_prime(z,alpha,prev)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
if alpha == 0
    y = 1/sqrt(z)*acosh((1+eta(0,z)-2*eta(0,z/4))/2/(eta(0,z/4)-1));
elseif alpha == 1
    y = fsolve(@(th) H(z,th),prev,opt);
else
    y = fsolve(@(th) G(z,alpha^2*z,th^2)-G(z,4*z,th^2),prev,opt);
end
end

function y=G(a,b,t2)
y=(eta(0,a/4)-eta(0,b/4))/(eta(-1,a*t2)-eta(-1,b*t2));
end
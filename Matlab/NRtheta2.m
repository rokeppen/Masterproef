% bepaling van de theta's voor RKs2
%   @param z: waarde voor mu1*h
%   @param alpha: verhouding z2/z1
function y = NRtheta2(z,alpha)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
if alpha == 1
    y = fzero(@(th) (eta(1,z/4)*eta(-1,z*th^2)/4-eta(0,z/4)*eta(0,z*th^2)*th^2)/2/eta(-1,z*th^2)^2, sqrt(3)/6);
else
    y = fsolve(@(th) eta(0,z/4)/eta(-1,z*th^2)-eta(0,alpha^2*z/4)/eta(-1,alpha^2*z*th^2), sqrt(3)/6, opt);
end
end
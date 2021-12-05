function y = NRtheta3(z,alpha,prev)
global opt;
if alpha == 0
    y = 1/sqrt(z)*acosh((1+eta(0,z)-2*eta(0,z/4))/2/(eta(0,z/4)-1));
elseif alpha == 1
    y = sqrt(15)/10;
else
    y = fsolve(@(th) H(z,alpha*z,th),prev,opt);
end
end

function y = H(z,z2,theta)
y=4*z2*(sinh(z)*(cosh(z*theta)^2 - cosh(z2*theta)/2 - 1/2)*sinh(z/2) - ((cosh(z) - 1)*(cosh(z) + 1)*(-cosh(z2*theta) + cosh(z*theta)))/4)*cosh(z2/2) - 4*sinh(z/2)*(cosh(z*theta) + 1/2)*sinh(z2)*z*cosh(z/2)*(cosh(z*theta) - 1);
end
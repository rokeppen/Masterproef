function y = NRtheta4(z,a1,a2,prev)
global opt;
fun = @(t) [b1(z,t(1),t(2)) - b1p(z,a1*z,t(1),t(2)); b1(z,t(1),t(2)) - b1p(z,a2*z,t(1),t(2))]; 
y = fsolve(fun,prev,opt);
end

function y = b1(z,th1,th2)
y=sinh(z/2)*(-cosh(2*z*th2)+cosh(z*th2)*cosh(z/2))/(z*(2*cosh(z*th1)*cosh(z*th2)*cosh(z/2)-cosh(z*th2)*cosh(2*z*th1)-cosh(z*th1)*cosh(2*z*th2)));
end

function y = b1p(z,z2,th1,th2)
y=-(sinh(z2/2)*cosh(z*th2)*z-sinh(z/2)*cosh(z2*th2)*z2)/(z*z2*(cosh(z*th1)*cosh(z2*th2)-cosh(z*th2)*cosh(z2*th1)));
end
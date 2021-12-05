function y = RKs3_fixed(xspan,h,alpha,y0,f,zfun)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
xn = xspan(1);
y = y0.';
sol = ones(1,4*length(y0))*0.1;
z = zfun(y(:,end),h);
th = NRtheta3(z^2,alpha^2*z^2);
[A,b,c] = calculate(z,th);
while xn < xspan(2)-h/2
    xn = xn+h;
    sol = fsolve(@(Y) RKs3_system(Y,A,b,c,h,f,y(:,end),xn),sol,opt);
    y = [y,sol(4:4:length(y0)*4).'];
end
end

function F = RKs3_system(Y,A,b,c,h,f,yn,xn)
Y1 = [xn-h+c(1)*h,Y(1:4:length(Y))];
Y2 = [xn-h+c(2)*h,Y(2:4:length(Y))];
Y3 = [xn-h+c(3)*h,Y(3:4:length(Y))];
fy = [f(Y1); f(Y2); f(Y3)];
F = zeros(1,length(Y));
for i=1:4:length(Y)
    j = (i-1)/4+1;
    F(i) = yn(j)+h*(A(1,1)*fy(1,j)+A(1,2)*fy(2,j)+A(1,3)*fy(3,j))-Y(i);
    F(i+1) = yn(j)+h*(A(2,1)*fy(1,j)+A(2,2)*fy(2,j)+A(2,3)*fy(3,j))-Y(i+1);
    F(i+2) = yn(j)+h*(A(3,1)*fy(1,j)+A(3,2)*fy(2,j)+A(3,3)*fy(3,j))-Y(i+2);
    F(i+3) = yn(j)+h*(b(1)*fy(1,j)+b(2)*fy(2,j)+b(1)*fy(3,j))-Y(i+3);
end
end

function [A,b,c] = calculate(z,th)
c = [1/2-th, 1/2, 1/2+th];
b1 = (eta(0,z^2)-eta(0,z^2/4))/2/(eta(-1,4*z^2*th^2)-eta(-1,z^2*th^2));
b2 = (eta(-1,4*z^2*th^2)*eta(0,z^2/4)-eta(0,z^2)*eta(-1,z^2*th^2))/(eta(-1,4*z^2*th^2)-eta(-1,z^2*th^2));
a2 = (cosh(2*z*th)-cosh(z*th)*cosh(z/2))/z/sinh(z*th);
a3 = (-cosh(z*th)+cosh(z/2))/z/sinh(z*th);
a4 = -b1*a2/b2;
A = [b1/2, b2/2-a2, b1/2-a3; b1/2-a4, b2/2, b1/2+a4; b1/2+a3, b2/2+a2, b1/2];
b = [b1,b2,b1];
end
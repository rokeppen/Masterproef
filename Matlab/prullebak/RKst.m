function y = RKst(xspan,h,alpha,y0,f,zfun)
xn = xspan(1);
y = y0.';
sol = ones(1,4*length(y0))*0.1;
while xn < xspan(2)-h/2
    z = feval(zfun,y(:,end),h);
    th = NRthetat(z^2,alpha,sqrt(15)/10);
    xn = xn+h;
    sol = fsolve(@(Y) RKs3_system(Y,z,h,f,y(:,end),th,xn),sol,optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16));
    subsol = [];
    for i = 1:length(y0)
        subsol = [subsol;sol(4*i)];
    end
    y = [y,subsol];
end
end

function z = consz(y,h)
z = 1i*h;
end

function z = conszr(y,h)
z = h;
end

function z = testz2(y,h)
global e
z = 1i*h*e;
end

function F = RKs3_system(Y,z,h,f,yn,th,xn)
c1 = 1/2-th;
c2 = 1/2;
c3 = 1/2+th;
b1 = (2*sinh(z/2)/z-1)/(2*cosh(th*z)-2);
b2 = (2*exp(th*z) - 2*exp(z*(th + 1)) + z*exp(z/2) + z*exp((z*(4*th + 1))/2))/(z*(exp(z/2) - 2*exp((z*(2*th + 1))/2) + exp((z*(4*th + 1))/2)));
a2 = 0;
a3 = (exp(2*th*z) - 1)/(z*(exp(2*th*z) + 1));
a4 = 0;
g1 = (2*exp(z/2) + 2*exp((z*(8*th + 1))/2))/(exp(z*(th + 1)) + exp(z*(3*th + 1)) + exp(th*z) + exp(3*th*z));
g2 = (2*exp(z/2))/(exp(z) + 1);
a11 = g1*b1/2;
a12 = g1*b2/2-a2;
a13 = g1*b1/2-a3;
a21 = g2*b1/2-a4;
a22 = g2*b2/2;
a23 = g2*b1/2+a4;
a31 = g1*b1/2+a3;
a32 = g1*b2/2+a2;
a33 = g1*b1/2;

Y1 = [xn-h+c1*h];
Y2 = [xn-h+c2*h];
Y3 = [xn-h+c3*h];
for i=1:4:length(Y)
    Y1 = [Y1,Y(i)];
    Y2 = [Y2,Y(i+1)];
    Y3 = [Y3,Y(i+2)];
end
fy1 = f(Y1);
fy2 = f(Y2);
fy3 = f(Y3);
for i=1:4:length(Y)
    j = (i-1)/4+1;
    F(i) = yn(j)+h*a11*fy1(j)+h*a12*fy2(j)+h*a13*fy3(j)-Y(i);
    F(i+1) = yn(j)+h*a21*fy1(j)+h*a22*fy2(j)+h*a23*fy3(j)-Y(i+1);
    F(i+2) = yn(j)+h*a31*fy1(j)+h*a32*fy2(j)+h*a33*fy3(j)-Y(i+2);
    F(i+3) = yn(j)+h*b1*fy1(j)+h*b2*fy2(j)+h*b1*fy3(j)-Y(i+3);
end
end

function y = NRthetat(z,alpha,prev)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
y = fsolve(@(th) Gt(z,th)-Gt(alpha^2*z,th),prev,opt);
end

function y = Gt(z,th)
y=(eta(0,z/4)-1)/2/(eta(-1,th^2*z)-1);
end
function y = RKs2_calvo(xspan,h,y0,f)
xn = xspan(1);
y = y0.';
sol = ones(1,3*length(y0))*0.1;
while xn < xspan(2)-h/2
    z = h;
    xn = xn+h;
    sol = fsolve(@(Y) RKs2_system(Y,z,h,f,y(:,end),xn),sol,optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16));
    subsol = [];
    for i = 1:length(y0)
        subsol = [subsol;sol(3*i)];
    end
    y = [y,subsol];
end
end

function F = RKs2_system(Y,z,h,f,yn,xn)
th = acos((cos(z/2)+sqrt(8+cos(z/2)^2))/4)/z;
c1 = 1/2-th;
c2 = 1/2+th;
g = 1;
b = sin(z/2)/z/cos(z*th);
a11 = (cos(2*z*th)-cos(z*th+z/2))/z/sin(2*z*th);
a22 = a11;
a12 = (cos(z*(th-1/2))-1)/z/sin(2*z*th);
a21 = a12;

Y1 = [xn-h+c1*h];
Y2 = [xn-h+c2*h];
for i=1:3:length(Y)
    Y1 = [Y1,Y(i)];
    Y2 = [Y2,Y(i+1)];
end
fy1 = f(Y1);
fy2 = f(Y2);
for i=1:3:length(Y)
    j = (i-1)/3+1;
    F(i) = g*yn(j)+h*a11*fy1(j)+h*a12*fy2(j)-Y(i);
    F(i+1) = g*yn(j)+h*a21*fy1(j)+h*a22*fy2(j)-Y(i+1);
    F(i+2) = yn(j)+h*b*fy1(j)+h*b*fy2(j)-Y(i+2);
end
end
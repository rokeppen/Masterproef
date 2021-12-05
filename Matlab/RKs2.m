% tweetraps EFRK-methode gefit op z1=zfun(h,y) en z2=alpha*z1
%   @param xspan: het interval waarover geïntegreerd wordt
%   @param h: de te gebruiken stapgrootte
%   @param alpha: de verhouding z2/z1
%   @param y0: de beginvoorwaarde
%   @param f: het probleemstelsel
%   @param zfun: de updatefunctie voor z1
function y = RKs2(xspan,h,alpha,y0,f,zfun)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
xn = xspan(1);
y = y0.';
sol = ones(1,3*length(y0))*0.1;
while xn < xspan(2)-h/2
    z = zfun(y(:,end),h);
    th = NRtheta2(z^2,alpha);
    xn = xn+h;
    sol = fsolve(@(Y) RKs2_system(Y,z,h,f,y(:,end),th,xn),sol,opt);
    y = [y,sol(3:3:length(y0)*3).'];
end
end

% impliciet stelsel ter bepaling van de volgende waarde van y
function F = RKs2_system(Y,z,h,f,yn,th,xn)
c = [1/2-th, 1/2+th];
g = cosh(2*z*th)/cosh(z/2)/cosh(z*th);
b = eta(0,z^2/4)/2/eta(-1,z^2*th^2);
l = -sinh(z*th)/cosh(z*th)/z;
A = [g*b/2, g*b/2+l; g*b/2-l, g*b/2];

Y1 = [xn-h+c(1)*h,Y(1:3:length(Y))];
Y2 = [xn-h+c(2)*h,Y(2:3:length(Y))];
fy = [f(Y1); f(Y2)];
for i=1:3:length(Y)
    j = (i-1)/3+1;
    F(i) = g*yn(j)+h*(A(1,1)*fy(1,j)+A(1,2)*fy(2,j))-Y(i);
    F(i+1) = g*yn(j)+h*(A(2,1)*fy(1,j)+A(2,2)*fy(2,j))-Y(i+1);
    F(i+2) = yn(j)+h*b*(fy(1,j)+fy(2,j))-Y(i+2);
end
end
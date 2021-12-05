function stability
omega = 1;
x = linspace(-10,5,1000);
y = linspace(-5,5,1000);
[X,Y] = meshgrid(x,y);
Z = abs(R(X+1i*Y,omega));
ZZ = abs(R(X+1i*Y,omega))-abs(exp(X+1i*Y));

figure
contour(X,Y,Z,[1,1],'Fill','on')

figure
contour(X,Y,ZZ,[0,0],'Fill','on')

%scatter(x,y,'red','filled');
%range = linspace(0,-20,101);
%for i=1:101
%    t = range(i);
%    tt(i)=R(t,-1);
%end
%figure;
%plot(range,tt);
end

function y = R(z,zf)
alpha = 2;
th = NRtheta(zf^2,alpha,sqrt(3)/6);
g = cosh(2*zf*th)/cosh(zf/2)/cosh(zf*th);
b = sinh(zf/2)/zf/cosh(zf*th);
l = -sinh(zf*th)/cosh(zf*th)/zf;
y = (1+g*b*z+l^2*z.^2)./(1-g*b*z+l^2*z.^2);
end


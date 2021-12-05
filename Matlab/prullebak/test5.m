function test5
h = 0.2;
z = 1i*h;
alpha = 2;
th = 0.288:0.00001:0.289;
y = sinh(z/2)./cosh(z*th)./z-sinh(alpha*z/2)./cosh(alpha*z*th)./alpha./z;
plot(th,y);

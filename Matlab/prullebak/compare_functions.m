function compare_functions
th = sqrt(15)/10;
zrange = linspace(0,5,501);
for i=1:501
    z = zrange(i)*1i;
    %y(i) = (sinh(z)-2*sinh(z/2))/2/z/(cosh(2*th*z)-cosh(th*z));
    k(i) = (eta(0,z^2)-eta(0,z^2/4))/(eta(-1,4*z^2*th^2)-eta(-1,z^2*th^2))/2;
    %y(i) = -(exp(2*z*(th + 1)) + 2*exp((z*(4*th + 1))/2) - 2*exp((z*(4*th + 3))/2) - exp(2*th*z))/(2*z*(exp(z*(th + 1)) + exp(z*(3*th + 1)) - exp(z*(4*th + 1)) - exp(z)));
    y(i)=-(exp(2*z*(th + 1)) + 2*exp((z*(4*th + 1))/2) - 2*exp((z*(4*th + 3))/2) - exp(2*th*z))/(2*z*(exp(z*(th + 1)) + exp(z*(3*th + 1)) - exp(z*(4*th + 1)) - exp(z)));
    k(i)-y(i)
end
figure;
plot(zrange,y,zrange,k)
legend('lang','etas');
end

%2*exp((z*(4*th + 3))/2) - 2*exp((z*(4*th + 1))/2) - exp(2*z*(th + 1)) + exp(2*th*z)
%2*z*(exp(z*(th + 1)) + exp(z*(3*th + 1)) - exp(z*(4*th + 1)) - exp(z))
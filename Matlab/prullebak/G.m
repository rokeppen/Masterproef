function y = G(z,th)
y=(eta(0,z)-eta(0,z/4))/(eta(-1,4*z*th^2)-eta(-1,z*th^2))/2;
%z=sqrt(z);
%y=-(exp(2*z*(th + 1)) + 2*exp((z*(4*th + 1))/2) - 2*exp((z*(4*th + 3))/2) - exp(2*th*z))/(2*z*(exp(z*(th + 1)) + exp(z*(3*th + 1)) - exp(z*(4*th + 1)) - exp(z)));
end
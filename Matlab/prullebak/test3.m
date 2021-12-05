function y = test3(z,th)
alpha = 2;
y = eta(0,z/4)/eta(-1,z*th^2)-eta(0,alpha^2*z/4)/eta(-1,alpha^2*z*th^2);
end
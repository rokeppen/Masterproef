function approxtest(diff,zmin,zmax)
in=501;
z1vals=linspace(zmin,zmax,in); %z1vals zijn kwadraten
for i=1:in
    z = z1vals(i);
    if z < 0
        alpha = imag((diff/i+sqrt(z1vals(i))))/sqrt(-z1vals(i));
    else
        alpha = abs((diff+sqrt(z1vals(i))))/sqrt(z1vals(i));
    end
    yapprox(i) = fzero(@(th) Fprime(z,th)+diff*Fdoubleprime(z,th)/2,sqrt(3)/6);
    diff*Fdoubleprime(z,yapprox(i))/2
    y(i) = NRtheta(z1vals(i),alpha,sqrt(3)/6);
    ydiff(i) = abs(y(i)-yapprox(i));
end
figure;
plot(z1vals,ydiff)
xlabel('Z_1')
ylabel('fout')
figTitleName = sprintf('Foutschatting bij Z_1 - Z_2 = %d',diff);
title(figTitleName)
end

function y = Fprime(z,th)
y = (eta(1,z/4)*eta(-1,z*th^2)/4-th^2*eta(0,z/4)*eta(0,z*th^2))/eta(-1,z*th^2)^2;
end

function y = Fdoubleprime(z,th)
y = ((1/4*(1/4*eta(2,z/4)*eta(-1,z*th^2)+eta(1,z/4)*th^2*eta(0,z*th^2))-th^2*(1/4*eta(1,z/4)*eta(0,z*th^2)+eta(0,z/4)*eta(1,z*th^2)*th^2))*eta(-1,z*th^2)^2-(eta(1,z/4)*eta(-1,z*th^2)/4-th^2*eta(0,z/4)*eta(0,z*th^2))*2*eta(-1,z*th^2)*eta(0,z*th^2)*th^2)/eta(-1,z^th^2)^4;
end
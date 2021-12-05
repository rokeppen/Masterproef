function test2
zrange = linspace(0,5,501);
for alpha = 0.1:0.1:3
    for i = 1:501
        z = zrange(i);
        k(i) = find_theta(z.^2,alpha^2*z.^2);
        y(i) = NRtheta3(z.^2,alpha,sqrt(15)/10);
        t(i) = k(i)-y(i);
    end
end
figure;
plot(zrange,t)

in = 51;
jn = 51;
z1vals=linspace(0,5,in);
z2vals=linspace(0,5,jn);
[X,Y] = meshgrid(z1vals,z2vals);
for i=1:in
    for j=1:jn
        k = find_theta(z1vals(i).^2,z2vals(j)^2*z1vals(i).^2);
        y = NRtheta3(z1vals(i).^2,z2vals(j),sqrt(15)/10);
        G(i,j) = k-y;
    end
end
figure;
surf(X,Y,G')
end
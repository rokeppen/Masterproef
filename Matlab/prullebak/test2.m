function test2(z)
alpha = 2;
thrange = linspace(0.1,0.9,501);
for i=1:501
    th = thrange(i);
    y(i) = eta(0,z/4)/eta(-1,z*th^2)-eta(0,alpha^2*z/4)/eta(-1,alpha^2*z*th^2);
end
figure;
plot(thrange,y)
end
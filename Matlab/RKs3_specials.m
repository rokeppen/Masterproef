% plot voor RKs3 het verloop van theta voor vier specifieke gevallen
function RKs3_specials
range = -5:0.01:5;
for i = 1:1001
    z1 = range(i);
    thz20(i) = NRtheta3(z1,0); % z2 = 0
    thz2z(i) = NRtheta3(z1,z1); % z2 = z1 
    thz23z(i) = NRtheta3(z1,3*z1); % z2 = 3*z1
    thz205z(i) = NRtheta3(z1,0.5*z1); % z2 = z1/2
end
figure;
plot(range,thz20);
hold on;
plot(range,thz2z);
plot(range,thz23z);
plot(range,thz205z);
hold off;
xlabel('Z_1')
ylabel('\theta')
legend('$z_2=0$','$z_2=z_1$','$z_2=3z_1$','$z_2=\frac{z_1}{2}$','Interpreter','latex','location','northwest');
end
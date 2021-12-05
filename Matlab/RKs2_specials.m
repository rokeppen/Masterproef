% plot voor RKs2 het verloop van theta voor drie specifieke gevallen
function RKs2_specials
range = -5:0.01:5;
for i = 1:1001
    z1 = range(i);
    thz20(i) = NRtheta2(z1,0); % z2 = 0
    thz2z(i) = NRtheta2(z1,1); % z2 = z1
    thz22z(i) = NRtheta2(z1,2); % z2 = 2*z1
end
figure;
plot(range,thz20);
hold on;
plot(range,thz2z);
plot(range,thz22z);
hold off;
xlabel('Z_1')
ylabel('\theta')
legend('$z_2=0$','$z_2=z_1$','$z_2=2z_1$','Interpreter','latex','location','northwest');
end
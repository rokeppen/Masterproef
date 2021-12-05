% plot de ligging van het gaussische geval voor theta en rondliggende
% waarden
function RKs3_cons_theta
th = sqrt(15)/10;
th2 = th + 0.001;
th3 = th - 0.001;
prec = 101;
z1 = linspace(-5,5,prec);
z2 = linspace(-5,5,prec);
for i = 1:prec
    for j = 1:prec
        Z(i,j) = G(z1(i),z2(j),th^2)-G(z1(i),4*z1(i),th^2);
        Z2(i,j) = G(z1(i),z2(j),th2^2)-G(z1(i),4*z1(i),th2^2);
        Z3(i,j) = G(z1(i),z2(j),th3^2)-G(z1(i),4*z1(i),th3^2);
    end
end
figure;
contour(z1,z2,Z',[0,0],'b');
hold on
contour(z1,z2,Z2',[0,0],'r');
contour(z1,z2,Z3',[0,0],'g');
xlabel('Z_1');
ylabel('Z_2');
axis([-5 5 -5 5])
legend('$\theta=\frac{\sqrt{15}}{10}$','$\theta=\frac{\sqrt{15}}{10}+\delta$','$\theta=\frac{\sqrt{15}}{10}-\delta$','Interpreter','latex')
hold off;
end
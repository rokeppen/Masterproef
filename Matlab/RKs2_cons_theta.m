% plot de ligging van het gaussische geval voor theta en rondliggende
% waarden
function RKs2_cons_theta
th = sqrt(3)/6;
th2 = th + 0.001;
th3 = th - 0.001;
prec = 101;
z1 = linspace(-5,5,prec);
z2 = linspace(-5,5,prec);
for i = 1:prec
    for j = i:prec
        Z(i,j) = eta(0,z1(i)/4)/eta(-1,z1(i)*th^2)-eta(0,z2(j)/4)/eta(-1,z2(j)*th^2);
        Z2(i,j) = eta(0,z1(i)/4)/eta(-1,z1(i)*th2^2)-eta(0,z2(j)/4)/eta(-1,z2(j)*th2^2);
        Z3(i,j) = eta(0,z1(i)/4)/eta(-1,z1(i)*th3^2)-eta(0,z2(j)/4)/eta(-1,z2(j)*th3^2);
        % gebruik maken van de symmetrie is efficiënter
        if i ~= j
            Z(j,i) = Z(i,j);
            Z2(j,i) = Z2(i,j);
            Z3(j,i) = Z3(i,j);
        end
    end
end
figure;
contour(z1,z2,Z',[0,0],'b');
hold on
contour(z1,z2,Z2',[0,0],'r');
contour(z1,z2,Z3',[0,0],'g');
plot(z1,-z1,'black');
xlabel('Z_1');
ylabel('Z_2');
legend('$\theta=\frac{\sqrt{3}}{6}$','$\theta=\frac{\sqrt{3}}{6}+\delta$','$\theta=\frac{\sqrt{3}}{6}-\delta$','$Z_2=-Z_1$','Interpreter','latex')
hold off;
end

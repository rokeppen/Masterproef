function cons_theta_RKs3_bis
A1=[-0.9 0];
B1=[4.6 0];
A2=[0.2 1.1];
B2=[4.6 0];
A3=[-2 -1.1];
figure;
axis([-5 5 -5 5])
xlim = get(gca,'XLim');
m = (B1(2)-B1(1))/(A1(2)-A1(1));
n = B1(2)*m - A1(2);
y1 = m*xlim(1) + n;
y2 = m*xlim(2) + n;
hold on
line([xlim(1) xlim(2)],[y1 y2],'Color','blue')
line([xlim(1)+1.1 xlim(2)+1.1],[y1 y2],'Color','red')
line([xlim(1)-1.1 xlim(2)-1.1],[y1 y2],'Color','green')
xlabel('Z_1');
ylabel('Z_2');
legend('$\theta=\frac{\sqrt{15}}{10}$','$\theta=\frac{\sqrt{15}}{10}+\delta$','$\theta=\frac{\sqrt{15}}{10}-\delta$','Interpreter','latex')
hold off;
end
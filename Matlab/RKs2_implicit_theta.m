% contourplot van de vergelijking ter bepaling van theta bij RKs2, met een
% doorsnede bij het geval z1
%   @param alpha: de verhouding z2/z1
%   @param z1: de waarde voor z1 ter bepaling van de doosnede
function RKs2_implicit_theta(alpha,z1)
prec = 501;
z = linspace(-15,0,prec);
th = linspace(0,1,prec);
for i = 1:prec
    for j = 1:prec
        Z(i,j) = real(eta(0,z(i)/4)./eta(-1,z(i).*th(j).^2)-eta(0,alpha^2*z(i)/4)./eta(-1,alpha^2*z(i).*th(j).^2));
    end
end
for j = 1:prec
    ZC(j) = real(eta(0,z1/4)./eta(-1,z1*th(j).^2)-eta(0,alpha^2*z1/4)./eta(-1,alpha^2*z1*th(j).^2));
end
figure;
contourf(z,th,Z');
hold on;
hRed = plot(NaN, '-g');
hBlue = plot(NaN, '-b');
xline(z1,'r','LineWidth',1.5);
xlabel('Z_1');
ylabel('\theta');
colormap(winter);
legend([hRed hBlue], '$F(Z_1)>F(Z_2)$', '$F(Z_1)<F(Z_2)$','Interpreter','latex','Location','northwest');
hold off;
figure;
plot(th,ZC,'r')
xlabel('\theta');
axis([0 1 -10 10])
end
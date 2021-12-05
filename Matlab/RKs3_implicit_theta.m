% contourplot van de vergelijking ter bepaling van theta bij RKs3, met een
% doorsnede bij het geval z1
%   @param alpha: de verhouding z2/z1
%   @param z1: de waarde voor z1 ter bepaling van de doosnede
function RKs3_implicit_theta(alpha,z1)
prec = 501;
z = linspace(-30,0,prec);
th = linspace(0,1,prec);
for i = 1:prec
    for j = 1:prec
        Z(i,j) = G(z(i),alpha^2*z(i),th(j)^2)-G(z(i),4*z(i),th(j)^2);
    end
end
for j = 1:prec
    ZC(j) = G(z1,alpha^2*z1,th(j)^2)-G(z1,4*z1,th(j)^2);
end
figure;
contourf(z,th,Z');
colormap(winter);
hold on;
hRed = plot(NaN, 'Color', [58/256 166/256 114/256]);
hBlue = plot(NaN, '-g');
xline(z1,'r','LineWidth',1.5);
xlabel('Z_1');
ylabel('\theta');
legend([hRed hBlue], '$G(z_1,z_2)>G(z_1,2z_1)$', '$G(z_1,z_2)<G(z_1,2z_1)$','Interpreter','latex','Location','northwest');
hold off;
figure;
plot(th,ZC,'r')
xlabel('\theta');
axis([0 1 -10 10])
end
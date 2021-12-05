% plot het verloop van de theta-coëfficiënten van de RKs3-methode in een
% vierkant interval op de fitting frequenties Z_i
%   @param zmin: benedengrens voor de Z_i
%   @param zmax: bovengrens voor de Z_i
function Z = RKs3_bepaal_thetas(zmin,zmax)
n = 51;
z1vals = linspace(zmin,zmax,n);
z2vals = linspace(zmin,zmax,n);
[X,Y] = meshgrid(z1vals,z2vals);
for i = 1:n
    for j = 1:n
        Z(i,j) = NRtheta3(z1vals(i),z2vals(j));
    end
end
figure
surf(X,Y,Z')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\theta')
axis([z1vals(1) z1vals(n) z2vals(1) z2vals(n) 0.38 0.395])
end
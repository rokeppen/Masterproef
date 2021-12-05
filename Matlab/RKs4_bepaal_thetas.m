% plot het verloop van de theta-coëfficiënten van de RKs4-methode in een
% vierkant interval op de fitting frequenties Z_i
%   @param zmin: benedengrens voor de Z_i
%   @param zmax: bovengrens voor de Z_i
function [Z1,Z2] = RKs4_bepaal_thetas(zmin,zmax)
n = 51;
z1vals = linspace(zmin,zmax,n);
z2vals = linspace(zmin,zmax,n);
[X,Y] = meshgrid(z1vals,z2vals);
for i = 1:n
    for j = i:n
        p = NRtheta4(z1vals(i),z2vals(j));
        Z1(i,j) = p(1);
        Z2(i,j) = p(2);
        if i ~= j
            Z1(j,i) = Z1(i,j);
            Z2(j,i) = Z2(i,j);
        end
    end
end
figure;
surf(X,Y,Z1')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\theta_1')

figure;
surf(X,Y,Z2')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\theta_2')
end
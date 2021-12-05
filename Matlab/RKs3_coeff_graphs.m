% plot het verloop van de coëfficiënten van de RKs3-methode in een
% vierkant interval op de fitting frequenties Z_i
%   @param zmin: benedengrens voor de Z_i
%   @param zmax: bovengrens voor de Z_i
function RKs3_coeff_graphs(zmin,zmax)
Z = RKs3_bepaal_thetas(zmin,zmax);
n = 51;
z1vals = linspace(zmin,zmax,n);
z2vals = linspace(zmin,zmax,n);
[X,Y] = meshgrid(z1vals,z2vals);
for i=1:n
    for j=1:n
        B1(i,j) = b1(z1vals(i),Z(i,j));
        B2(i,j) = b2(z1vals(i),Z(i,j));
        A2(i,j) = alpha2(z1vals(i),Z(i,j));
        A3(i,j) = alpha3(z1vals(i),Z(i,j));
    end
end
figure;
surf(X,Y,B1')
xlabel('Z_1')
ylabel('Z_2')
zlabel('b_1')

figure;
surf(X,Y,B2')
xlabel('Z_1')
ylabel('Z_2')
zlabel('b_2')

figure;
surf(X,Y,A2')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\alpha_2')

figure;
surf(X,Y,A3')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\alpha_3')
end

% b_1-coëfficiënt van RKs3
function y = b1(z,th)
if z == 0
    y = 1/(24*th^2);
else
    y = (eta(0,z)-eta(0,z/4))/2/(eta(-1,4*z*th^2)-eta(-1,z*th^2));
end
end

% b_2-coëfficiënt van RKs3
function y = b2(z,th)
if z == 0
    y = 1 - 1/(12*th^2);
else
    y = (eta(-1,4*z*th^2)*eta(0,z/4)-eta(0,z)*eta(-1,z*th^2))/(eta(-1,4*z*th^2)-eta(-1,z*th^2));
end
end

% alpha_2-coëfficiënt van RKs3
function y = alpha2(z,th)
if z == 0
    y = (3*th)/2 - 1/(8*th);
else
    y = (cosh(2*z*th)-cosh(z*th)*cosh(z/2))/z/sinh(z*th);
end
end

% alpha_3-coëfficiënt van RKs3
function y = alpha3(z,th)
if z == 0
    y = -th/2 + 1/(8*th);
else
    y = (-cosh(z*th)+cosh(z/2))/z/sinh(z*th);
end
end
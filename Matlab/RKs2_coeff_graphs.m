% plot het verloop van de coëfficiënten van de RKs2-methode in een
% vierkant interval op de fitting frequenties Z_i
%   @param zmin: benedengrens voor de Z_i
%   @param zmax: bovengrens voor de Z_i
function RKs2_coeff_graphs(zmin,zmax)
Z = RKs2_bepaal_thetas(zmin,zmax);
n = 51;
z1vals = linspace(zmin,zmax,n);
z2vals = linspace(zmin,zmax,n);
[X,Y] = meshgrid(z1vals,z2vals);
for i=1:n
    for j=1:n
        G(i,j)=gamma(z1vals(i),Z(i,j));
        B(i,j)=b(z1vals(i),Z(i,j));
        L(i,j)=lambda(z1vals(i),Z(i,j));
    end
end
figure;
surf(X,Y,G')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\gamma')

figure;
surf(X,Y,B')
xlabel('Z_1')
ylabel('Z_2')
zlabel('b')

figure;
surf(X,Y,L')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\lambda')
end

% gamma-coëfficiënt van RKs2
function y = gamma(z,th)
z = sqrt(z);
y = cosh(2*z*th)/cosh(z/2)/cosh(z*th);
end

% lambda-coëfficiënt van RKs2
function y = lambda(z,th)
if z == 0
    y = -th;
else
    z = sqrt(z);
    y = -sinh(z*th)/z/cosh(z*th);
end
end

% b-coëfficiënt van RKs2
function y = b(z,th)
if z == 0
    y = 1/2;
else
    z = sqrt(z);
    y = sinh(z/2)/z/cosh(z*th);
end
end
% plot het verloop van de b-coëfficiënten van de RKs4-methode in een
% vierkant interval op de fitting frequenties Z_i
%   @param zmin: benedengrens voor de Z_i
%   @param zmax: bovengrens voor de Z_i
function RKs4_coeff_graphs(zmin,zmax)
[Z1,Z2] = RKs4_bepaal_thetas(zmin,zmax);
n = 51;
z1vals = linspace(zmin,zmax,n);
z2vals = linspace(zmin,zmax,n);
[X,Y] = meshgrid(z1vals,z2vals);
for i=1:n
    for j=1:n
        B1(i,j) = b1(z1vals(i),z2vals(j),Z1(i,j),Z2(i,j));
        B2(i,j) = b2(z1vals(i),z2vals(j),Z1(i,j),Z2(i,j));
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
end

% b_1-coëfficiënt van RKs4
function y = b1(z,z2,t1,t2)
if z == 0 && z2 == 0
    y = 2*(1/8-sqrt(30)/144);
elseif z == 0 || z2 == 0
    y = (-12*t2^2 + 1)/(24*t1^2 - 24*t2^2);
elseif abs(z) == abs(z2)
    y = ((-2*t2*(z*cosh(z*t1)*t1-sinh(z*t1))*cosh(z*t2)^2+2*sinh(z*t2)*(sinh(z*t1)*t2^2*z-t1*cosh(z*t1))*cosh(z*t2)+2*cosh(z*t1)*t1*t2*z)*sinh(z/2)+z*cosh(z/2)*cosh(z*t2)*(t1*cosh(z*t1)*sinh(z*t2)-t2*sinh(z*t1)*cosh(z*t2)))/(2*(-2*t1*(cosh(z*t1)^2-1/2)*t2*cosh(z*t2)^2+cosh(z*t1)*sinh(z*t1)*sinh(z*t2)*(t1^2+t2^2)*cosh(z*t2)+cosh(z*t1)^2*t1*t2)*z^2);
else
    t1 = t1^2;
    t2 = t2^2;
    y = (eta(-1,z2*t2)*eta(0,z/4)-eta(-1,z*t2)*eta(0,z2/4))/2/(eta(-1,z2*t2)*eta(-1,z*t1)-eta(-1,z2*t1)*eta(-1,z*t2));
end
end

% b_2-coëfficiënt van RKs4
function y = b2(z,z2,t1,t2)
if z == 0 && z2 == 0
    y = 2*(1/8+sqrt(30)/144);
else
    y = b1(z,z2,t2,t1);
end
end
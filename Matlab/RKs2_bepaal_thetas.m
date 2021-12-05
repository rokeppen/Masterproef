% plot het verloop van de theta-coëfficiënten van de RKs2-methode in een
% vierkant interval op de fitting frequenties Z_i
%   @param zmin: benedengrens voor de Z_i
%   @param zmax: bovengrens voor de Z_i
function Z = RKs2_bepaal_thetas(zmin,zmax)
n=51;
z1vals=linspace(zmin,zmax,n);
z2vals=linspace(zmin,zmax,n);
[X,Y] = meshgrid(z1vals,z2vals);
for i=1:n
    for j=i:n
        if z1vals(i) ~= 0
            Z(i,j)=NRtheta2(z1vals(i),sqrt(z2vals(j)/z1vals(i)));
        elseif z2vals(j) == 0
            Z(i,j)=sqrt(3)/6;
        else
            Z(i,j)=NRtheta2(z2vals(j),0);
        end
        if i ~= j
            Z(j,i)=Z(i,j);
        end
    end
end
figure
surf(X,Y,Z')
xlabel('Z_1')
ylabel('Z_2')
zlabel('\theta')
axis([z1vals(1) z1vals(n) z2vals(1) z2vals(n) 0.28 0.30])
end
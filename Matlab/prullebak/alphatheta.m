function alphatheta(zmin,zmax,coeff,a)
in=51;
z1vals=linspace(zmin,zmax,in);
prev=sqrt(3)/6;
for i=in:-1:1
    %if i ~= in
    %    prev = Z4(i+1);
    %end
    %Z1(i)=NRtheta(z1vals(i),0,prev);
    %Z2(i)=NRtheta(z1vals(i),1/2,prev);
    %Z3(i)=NRtheta(z1vals(i),1,prev);
    Z4(i)=NRtheta(z1vals(i),a,prev);
    %Z5(i)=NRtheta(z1vals(i),3,prev);
end
figure;
%plot(z1vals, Z1, z1vals, Z2, z1vals, Z3, z1vals, Z4);
%legend('\alpha = 0','\alpha = 1/2','\alpha = 1','\alpha = 2','\alpha = 3', 'Location', 'southeast');
plot(z1vals,Z4);
title('\theta voor z_2 = \alpha z_1')
xlabel('Z = z^2')
ylabel('\theta')

if coeff
    for i=1:in
        %G1(i)=gamma(z1vals(i),Z1(i));
        %G2(i)=gamma(z1vals(i),Z2(i));
        %G3(i)=gamma(z1vals(i),Z3(i));
        G4(i)=gamma(z1vals(i),Z4(i));
        %G5(i)=gamma(z1vals(i),Z5(i));
    end
    figure;
    %plot(z1vals, G1, z1vals, G2, z1vals, G3, z1vals, G4);
    %legend('\alpha = 0','\alpha = 1/2','\alpha = 1','\alpha = 2','\alpha = 3', 'Location', 'southeast');
    plot(z1vals,G4);
    title('\gamma voor z_2 = \alphaz_1')
    xlabel('Z = z^2')
    ylabel('\gamma')

    for i=1:in
        %B1(i)=b1(z1vals(i),Z1(i));
        %B2(i)=b1(z1vals(i),Z2(i));
        %B3(i)=b1(z1vals(i),Z3(i));
        B4(i)=b1(z1vals(i),Z4(i));
        %B5(i)=b1(z1vals(i),Z5(i));
    end
    figure;
    %plot(z1vals, B1, z1vals, B2, z1vals, B3, z1vals, B4);
    %legend('\alpha = 0','\alpha = 1/2','\alpha = 1','\alpha = 2','\alpha = 3', 'Location', 'southeast');
    plot(z1vals,B4);
    title('b voor z_2 = \alphaz_1')
    xlabel('Z = z^2')
    ylabel('b')

    for i=1:in
        %L1(i)=lambda(z1vals(i),Z1(i));
        %L2(i)=lambda(z1vals(i),Z2(i));
        %L3(i)=lambda(z1vals(i),Z3(i));
        L4(i)=lambda(z1vals(i),Z4(i));
        %L5(i)=lambda(z1vals(i),Z5(i));
    end
    figure;
    %plot(z1vals, L1, z1vals, L2, z1vals, L3, z1vals, L4);
    %legend('\alpha = 0','\alpha = 1/2','\alpha = 1','\alpha = 2','\alpha = 3', 'Location', 'southeast');
    plot(z1vals,L4);
    title('\lambda voor z_2 = \alphaz_1')
    xlabel('Z = z^2')
    ylabel('\lambda')
end
end

function y = gamma(z,th)
z = sqrt(z);
y = cosh(2*z*th)/cosh(z/2)/cosh(z*th);
end

function y = lambda(z,th)
if z == 0
    y = -th;
else
    z = sqrt(z);
    y = -sinh(z*th)/z/cosh(z*th);
end
end

function y = b1(z,th)
if z == 0
    y = 1/2;
else
    z = sqrt(z);
    y = sinh(z/2)/z/cosh(z*th);
end
end
function stability2
x = linspace(-5,0,2);
y = linspace(0,5,2);
[X,Y] = meshgrid(x,y);
%Z = abs(1+(X+1i*Y)*sin(2)/2);
Z = abs(R(X+1i*Y))

figure
contour(X,Y,Z,[1,1],'Fill','on')
end

function y = R(z)
a = 0.5;
b = 1/12;
y = (1+a*z+b*z^2)./(1-a*z+b*z^2);
end

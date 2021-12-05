function test_find_theta;

%find_theta(1,3)
%find_theta2(-1,-9)
% z1=1;
% thetavals=[];
% for j=1:100;
%   thetavals=[thetavals,find_theta(z,z2vals(j))];
% end
% plot(z2vals,thetavals)

z1vals=linspace(-10,10,10);
z2vals=linspace(-8,8,10);
[X,Y] = meshgrid(z1vals,z2vals);
[m,n]=size(X)
[m,n]=size(Y)
for i=1:10
  for j=1:10;
  Z(i,j)=find_theta2(z1vals(i),z2vals(j));
end
end
Z
figure
surf(X,Y,Z')

find_theta(50,-50)

end

function th=find_theta(z,z2);
th=sqrt(15)/10;
G=inline('(sinh(a/2)/(a/2)-sinh(b/2)/(b/2))*cosh(c*theta)','a','b','c','theta');
H=inline('(sinh(a/2)/(a/2)-sinh(b/2)/(b/2))*sinh(c*theta)*c','a','b','c','theta');
del=1;
iter=0;
while abs(del)>10^(-5) 
   del=-(G(z,2*z,z2,th)+G(2*z,z2,z,th)+G(z2,z,2*z,th))/(H(z,2*z,z2,th)+H(2*z,z2,z,th)+H(z2,z,2*z,th));
   th=th+del;
   iter=iter+1;
   if iter>5 
    disp ('oeie') 
   end
end
end

function th=find_theta2(z,z2); % z'en zijn kwadraten !
th=sqrt(15)/10;
%G=inline('(sinh(a/2)/(a/2)-sinh(b/2)/(b/2))*cosh(c*theta)','a','b','c','theta');
%H=inline('(sinh(a/2)/(a/2)-sinh(b/2)/(b/2))*sinh(c*theta)*c','a','b','c','theta');
del=1;
iter=0;
while abs(del)>10^(-5)
   th2=th^2; 
   del=-(G(z,4*z,z2,th2)+G(4*z,z2,z,th2)+G(z2,z,4*z,th2))/(H(z,4*z,z2,th2)+H(4*z,z2,z,th2)+H(z2,z,4*z,th2));
   th=th+del;
   iter=iter+1;
   
end
if iter>5 
    fprintf ('oeie %d \n',iter) 
   end
end

function y=G(a,b,c,t2)
y=(eta(a/4)-eta(b/4))*xi(c*t2);
end

function y=H(a,b,c,t2) % t2 staat voor theta^2
if t2 > 0
  t=sqrt(t2);
else
  t=sqrt(-t2);
end;
y=(eta(a/4)-eta(b/4))*eta(c*t2)*c*t;
end
   
function y=xi(z)
if z>0
   y=cosh(sqrt(z));
elseif z<0
   y=cos(sqrt(-z));
else
   y=1;
end
end

function y=eta(z)
if z>0
   y=sinh(sqrt(z))/sqrt(z);
elseif z<0
   y=sin(sqrt(-z))/sqrt(-z);
else y=1;
end
end


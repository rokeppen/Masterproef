function probleem;

in=20;
jn=20;
z1vals=linspace(-20,20,in);
z2vals=linspace(-20,20,jn);
[X,Y] = meshgrid(z1vals,z2vals);
[m,n]=size(X);
[m,n]=size(Y);
for i=1:in
  for j=1:jn;
  Z(i,j)=find_theta(z1vals(i),z2vals(j));
end
end
figure
surf(X,Y,Z')
xlabel('z')
ylabel('z2')
axis([z1vals(1) z1vals(in) z2vals(1) z2vals(jn) 0 1])

z1vals(9),z2vals(20)
find_theta(z1vals(9),z2vals(20))

end

function th=find_theta(z,z2); % z'en zijn kwadraten !
d2 = 0.01;
if (abs(z)>d2) && (abs(z2)>d2) % geen reeksontwikkeling
    th=sqrt(15)/10;
    del=1;
    iter=0;
    if abs(z2-4*z)<0.0001
       z1=z2;
    else
       z1=z;
    end
    while abs(del)>10^(-5) && (iter<5)
       th2=th^2;
       del=-(G(z1,z2,th2)-G(z,4*z,th2))/(-G(z1,z2,th2)*th*H(z1,z2,th2)+G(z,4*z,th2)*th*H(z,4*z,th2));
       th=th+del;
       [z1,z2,th2,G(z1,z2,th2),G(z,4*z,th2),H(z1,z2,th2),H(z,4*z,th2),del,iter+1]
       iter=iter+1;
    end
    if iter>5 
        fprintf ('oeie %d \n',iter)
        th 
    end
else % geval van reeksontwikkelingen of z=0 of z2=0
    if (z2==0) && (abs(z)>d2) % expliciete berekening mogelijk
        th=1/z*arccosh((1+eta(z)-2*eta(z/4))/2/(eta(z/4)-1));
    elseif (z2<d2) && (z<d2)
        th=sqrt(15)/10* ...
            (1+(1/2100+(-131/105840000+(13487/4889808000000+ ...
            (-1175117/320380220160000000-505147/91537205760000000000*z2)*z2)*z2)*z2)*z2 ...
            +(1/420+(-17/21168000+(-1153/1955923200000+(-21517/16019011008000000+137477/4660075929600000000*z2)*z2)*z2)*z2...
            +(-17/784000+(-330733/1955923200000+(223250821/320380220160000000-43658273/25630417612800000000*z2)*z2)*z2...
            +(769/4346496000+(19182533/5339670336000000-8992463/1708694507520000000*z2)*z2...
            +(-280393/284782417920000-367176113/10252167045120000000*z2-7963/2531399270400000*z)*z)*z)*z)*z);        
    end    
end
end

function y=G(a,b,t2)
d=0.01;
if abs(a-b)<0.01 
   y=(xi(a/4)-eta(a/4))/a/t2/eta(t2*a);
else   
  y=(eta(a/4)-eta(b/4))/(xi(a*t2)-xi(b*t2));
end
end

function y=H(a,b,t2) % t2 staat voor theta^2
d=0.01;
if abs(a-b)<0.01 
   y=(eta(a*t2)+xi(a*t2))/t2/eta(t2*a);
else
   y=(a*eta(a*t2)-b*eta(b*t2))/(xi(a*t2)-xi(b*t2));
end;
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


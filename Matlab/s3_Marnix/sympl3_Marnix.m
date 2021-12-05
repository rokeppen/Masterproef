function [t,y]=sympl3_Marnix(f,tspan,y0,jac,h)
yn=y0;
tn=tspan(1);
tend=tspan(2);
t=tn;
y=yn;
while tn<tend
  yn=sympl3_Marnix_step(f,jac,yn,h);
  tn=tn+h;
  t=[t;tn];
  y=[y;yn];
  end
end

function y=sympl3_Marnix_step(f,jacob,yn,h)
global fac
if strcmp(func2str(f),'kepler')
  omega1=(yn(1)^2+yn(2)^2)^(-3/2);
  omega2=fac*(yn(1)^2+yn(2)^2)^(-3/2);
elseif strcmp(func2str(f),'pert_kepler')
  omega1=1;
  omega2=fac;
elseif strcmp(func2str(f),'euler_prob')
   omega1=2*pi/7.45056320933095;
   omega2=fac*2*pi/7.45056320933095;
elseif strcmp(func2str(f),'simple')
   omega1=1;
   omega2=fac;
end
z=sqrt(-1)*omega1*h;
z2=sqrt(-1)*omega2*h; 
if abs(z)>0.1
 theta=find_theta(z^2,z2^2);
 b1=-(sinh(z)-2*sinh(z/2))/2/z/(cosh(theta*z)-cosh(2*theta*z));
 b2=(cosh(theta*z)*sinh(z)-2*cosh(2*theta*z)*sinh(z/2))/z/(cosh(theta*z)-cosh(2*theta*z));
 g1=(1-cosh(z/2))/z/sinh(theta*z);
 g2=(cosh(theta*z)-cosh(z/2))/z/sinh(theta*z);
else
 sq=sqrt(15);
 theta=sq/10*(1+z^2/150-31*z^4/240000+89*z^6/144e6+45539*z^8/72576e7-3085681*z^10/145152e10);
 b1=5/18-1/270*z^2-23/432000*z^4+1433/226800000*z^6-555073/2612736000000*z^8+24846889/14370048000000000*z^10;
 b2=4/9+1/135*z^2+23/216000*z^4+37/7087500*z^6-216047/1306368000000*z^8+14276111/7185024000000000*z^10;
 g1=sq*(-1/12+13/14400*z^2+1/288000*z^4-28061/48384000000*z^6+1192963/87091200000000*z^8+1695787/23950080000000000*z^10);
 g2=sq*(-1/30+11/18000*z^2-11/1800000*z^4-6653/30240000000*z^6+107593/8709120000000*z^8-48160367/239500800000000000*z^10);
end    

% b1=5/18-1/270*z^2-23/432000*z^4+1433/226800000*z^6-555073/2612736000000*z^8+24846889/14370048000000000*z^10+73678387217/448345497600000000000*z^12-3694434776297/470762772480000000000000*z^14
% b2=4/9+1/135*z^2+23/216000*z^4+37/7087500*z^6-216047/1306368000000*z^8+14276111/7185024000000000*z^10+24421960183/224172748800000000000*z^12-1464866821153/235381386240000000000000*z^14
% g1=sq*(-1/12+13/14400*z^2+1/288000*z^4-28061/48384000000*z^6+1192963/87091200000000*z^8+1695787/23950080000000000*z^10-41892939139/2988969984000000000000*z^12+20415580357/61802702438400000000000*z^14+15815563300037/3104181190656000000000000000*z^16
% g2=sq*(-1/30+11/18000*z^2-11/1800000*z^4-6653/30240000000*z^6+107593/8709120000000*z^8-48160367/239500800000000000*z^10-3871172351/679311360000000000000*z^12+221449633549/570621542400000000000000*z^14-247569477511723/42682491371520000000000000000*z^16

a11=b1/2;
a12=b2/2/b1*(b1+g1);
a13=b1/2+g2;
a21=(b1-g1)/2;
a22=b2/2;
a23=(b1+g1)/2;
a31=b1/2-g2;
a32=b2/2/b1*(b1-g1);
a33=b1/2;
b3=b1;

Yn=[yn,yn,yn]';
[m,n]=size(yn);
del=1;
jac=feval(jacob,yn);
J=eye(3*n,3*n)-h*[a11*jac,a12*jac,a13*jac;a21*jac,a22*jac,a23*jac;a31*jac,a32*jac,a33*jac];
Y=[yn',yn',yn'];
while norm(del,2) > 10^(-12)
   Yv=[Y(:,1);Y(:,2);Y(:,3)];
   F=Yv-Yn-h*[a11*feval(f,Y(:,1))+a12*feval(f,Y(:,2))+a13*feval(f,Y(:,3));
             a21*feval(f,Y(:,1))+a22*feval(f,Y(:,2))+a23*feval(f,Y(:,3));
             a31*feval(f,Y(:,1))+a32*feval(f,Y(:,2))+a33*feval(f,Y(:,3))];
   del=J\(-F);
   Y=Y+[del(1:n),del(n+1:2*n),del(2*n+1:3*n)];        
end    
y=yn+h*((b1*feval(f,Y(:,1))+b2*feval(f,Y(:,2))+b3*feval(f,Y(:,3)))');
end

function th=find_theta(z,z2); % z'en zijn kwadraten !
d2 = 0.01;
if (z2==0) && (abs(z)>d2) % expliciete berekening mogelijk
        th=1/sqrt(z)*acosh((1+eta(z)-2*eta(z/4))/2/(eta(z/4)-1));
elseif (z==0) && (abs(z2)> d2)
    th=sqrt(15)/10;
    del=1;
    iter=0;
    while (abs(del)>10^(-10)) && (iter<10)
       th2=th^2;
       del=  -((xi(z2*th2)-1)/12/th2+1-eta(z2/4))/...
              (z2/12*eta(z2*th2)-(xi(z2*th2)-1)/6/th2)*th;
       th=th+del;
       iter=iter+1;
    end
elseif (abs(z)>d2) | (abs(z2)>d2) % geen reeksontwikkeling
    th=sqrt(15)/10;
    del=1;
    iter=0;
    if abs(z2-4*z)<0.0001
       z1=z2;
    else
       z1=z;
    end
    while (abs(del)>10^(-10)) && (iter<10)
       th2=th^2;
       del=-(G(z1,z2,th2)-G(z,4*z,th2))/(-G(z1,z2,th2)*th*H(z1,z2,th2)+G(z,4*z,th2)*th*H(z,4*z,th2));
       th=th+del;
       iter=iter+1;
    end
    if iter>5 
        fprintf ('oeie %d \n',iter)
        th 
    end
else % geval van reeksontwikkelingen 
    th=sqrt(15)/10* ...
            (1+(1/2100+(-131/105840000+(13487/4889808000000+ ...
            (-1175117/320380220160000000-505147/91537205760000000000*z2)*z2)*z2)*z2)*z2 ...
            +(1/420+(-17/21168000+(-1153/1955923200000+(-21517/16019011008000000+137477/4660075929600000000*z2)*z2)*z2)*z2...
            +(-17/784000+(-330733/1955923200000+(223250821/320380220160000000-43658273/25630417612800000000*z2)*z2)*z2...
            +(769/4346496000+(19182533/5339670336000000-8992463/1708694507520000000*z2)*z2...
            +(-280393/284782417920000-367176113/10252167045120000000*z2-7963/2531399270400000*z)*z)*z)*z)*z);        
       
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
end
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

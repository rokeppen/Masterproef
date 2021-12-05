function [t,y]=sympl3_Calvo_fixed(f,tspan,y0,jac,h)
yn=y0;
tn=tspan(1);
tend=tspan(2);
t=tn;
y=yn;
while tn<tend
  yn=sympl3_Calvo_fixed_step(f,jac,yn,h);
  tn=tn+h;
  t=[t;tn];
  y=[y;yn];
  end
end

function y=sympl3_Calvo_fixed_step(f,jacob,yn,h)

if strcmp(func2str(f),'kepler')
  omega=(yn(1)^2+yn(2)^2)^(-3/2);
elseif strcmp(func2str(f),'pert_kepler')
  omega=1;
elseif strcmp(func2str(f),'euler_prob')
   omega=2*pi/7.45056320933095;
end
z=sqrt(-1)*omega*h; 
if abs(z)>0.01
 theta=sqrt(15)/10;
 b1=-(sinh(z)-2*sinh(z/2))/2/z/(cosh(theta*z)-cosh(2*theta*z));
 b2=(cosh(theta*z)*sinh(z)-2*cosh(2*theta*z)*sinh(z/2))/z/(cosh(theta*z)-cosh(2*theta*z));
 g1=(1-cosh(z/2))/z/sinh(theta*z);
 g2=(cosh(theta*z)-cosh(z/2))/z/sinh(theta*z);
else
 sq=sqrt(15);
 b1=5/18+1/14400*z^4-191/87091200*z^6+623/8294400000*z^8-78713/30656102400000*z^10;
 b2=4/9-1/7200*z^4-241/43545600*z^6+217/4147200000*z^8-8147/3065610240000*z^10;
 g1=sq*(-1/12+1/2880*z^2-13/1728000*z^4+221/1935360000*z^6-6061/3483648000000*z^8+9733/367873228800000*z^10);
 g2=sq*(-1/30-1/3600*z^2+1/540000*z^4-17/604800000*z^6+23/54432000000*z^8-4603/718502400000000*z^10);
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
while max(norm(del,2)) > 10^(-15)
   Yv=[Y(:,1);Y(:,2);Y(:,3)];
   FY1=feval(f,Y(:,1));
   FY2=feval(f,Y(:,2));
   FY3=feval(f,Y(:,3));
   F=Yv-Yn-h*[a11*FY1+a12*FY2+a13*FY3;
             a21*FY1+a22*FY2+a23*FY3;
             a31*FY1+a32*FY2+a33*FY3]
   del=J\(-F);
   Y=Y+[del(1:n),del(n+1:2*n),del(2*n+1:3*n)];        
end
   
y=yn+h*((b1*feval(f,Y(:,1))+b2*feval(f,Y(:,2))+b3*feval(f,Y(:,3)))');
end
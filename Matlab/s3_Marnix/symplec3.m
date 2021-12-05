function fout=symplec3(h,methode);

tspan = [0, 1000];
e=0.001;
y0 = [1-e,0,0,sqrt((1+e)/(1-e))];

if methode==1
[t,y] = sympl3_klassiek(@kepler, tspan, y0 , @kepler_jac,h);
elseif methode==2
[t,y] = sympl3_Calvo(@kepler, tspan, y0 , @kepler_jac,h);
end
yend=y(end,1:2);
ysol=kepler_sol(tspan(2),e);
fout=norm(yend-ysol,2);

end

function [t,y]=sympl3_klassiek(f,tspan,y0,jac,h)
yn=y0;
tn=tspan(1);
tend=tspan(2);
t=tn;
y=yn;
while tn<tend
  yn=sympl3_klassiek_step(f,jac,yn,h);
  tn=tn+h;
  t=[t;tn];
  y=[y;yn];
  end
end

function y=sympl3_klassiek_step(f,jacob,yn,h)
sq=sqrt(15);
a11=5/36;
a12=2/9-sq/15;
a13=5/36-sq/30;
a21=5/36+sq/24;
a22=2/9;
a23=5/36-sq/24;
a31=5/36+sq/30;
a32=2/9+sq/15;
a33=5/36;
b1=5/18;
b2=4/9;
b3=5/18;

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
%-----------------------------------------------------------------------
function [t,y]=sympl3_Calvo(f,tspan,y0,jac,h)
yn=y0;
tn=tspan(1);
tend=tspan(2);
t=tn;
y=yn;
while tn<tend
  yn=sympl3_Calvo_step(f,jac,yn,h);
  tn=tn+h;
  t=[t;tn];
  y=[y;yn];
  end
end

function y=sympl3_Calvo_step(f,jacob,yn,h)
omega=(yn(1)^2+yn(2)^2)^(-3/2);
z=sqrt(-1)*omega*h;
if abs(z)>0.01
 beta=sqrt(5+2*cosh(z/2)+sqrt(15+8*cosh(z/2)+2*cosh(z)))/2/sqrt(3);
 theta=2*acosh(beta)/z;
 b1=-(sinh(z)-2*sinh(z/2))/2/z/(cosh(theta*z)-cosh(2*theta*z));
 b2=(cosh(theta*z)*sinh(z)-2*cosh(2*theta*z)*sinh(z/2))/z/(cosh(theta*z)-cosh(2*theta*z));
 g1=(1-cosh(z/2))/z/sinh(theta*z);
 g2=(cosh(theta*z)-cosh(z/2))/z/sinh(theta*z);
else
 sq=sqrt(15);
 theta=sq/10*(1+z^2/150-31*z^4/240000+89*z^6/144e6+45539*z^8/72576e7-3085681*z^10/145152e10);
 b1=5/18-z^2/270-23*z^4/432000+1433*z^6/2268e5-555073*z^8/2612736e6+24846889*z^10/14370048e9;
 b2=4/9+z^2/135+23*z^4/216000+37*z^6/7087500-216047*z^8/1306368e6+14276111*z^10/7185024e9;
 g1=sq*(-1/12+13*z^2/14400+z^4/288000-28061*z^6/48384e6+1192963*z^8/870912e8+1695787*z^10/2395008e10);
 g2=sq*(-1/30+11*z^2/18000-11*z^4/180000-6653/3024e7*z^6+107593*z^8/870912e7-48160367*z^10/2395008e11);
end    

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






function dydt = kepler(y)
noem=(y(1)^2+y(2)^2)^(3/2);
dydt = [y(3);
        y(4);
        -y(1)/noem;
        -y(2)/noem];
end
        
function dfdy = kepler_jac(y)    
noem3=(y(1)^2+y(2)^2)^(3/2);
noem5=(y(1)^2+y(2)^2)^(5/2);
dfdy= [0,0,1,0;
       0,0,0,1;
       -1/noem3+3*y(1)^2/noem5,3*y(1)*y(2)/noem5,0,0;
       3*y(1)*y(2)/noem5,-1/noem3+3*y(2)^2/noem5,0,0];
end

function sol=kepler_sol(t,e)
E=t;
del=1;
while abs(del)>10^(-12)
  del=-(t-E+e*sin(E))/(-1+e*cos(E));
  E=E+del;
end
sol=[cos(E)-e,sqrt(1-e^2)*sin(E)];
end


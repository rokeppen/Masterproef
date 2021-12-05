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
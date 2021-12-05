function sol=kepler_sol(t,e)
E=t;
del=1;
while abs(del)>10^(-10)
  del=-(t-E+e*sin(E))/(-1+e*cos(E));
  E=E+del;
end
c=cos(E);
s=sin(E);
noem=1/(1-e*c);
sq=sqrt(1-e^2);
sol=[c-e,sq*s,-s*noem,sq*noem*c];
end
function sol=pert_kepler_sol(t)

global ep;
a=1+ep;
c=cos(a*t);
s=sin(a*t);
sol=[c,s,-a*s,a*c];
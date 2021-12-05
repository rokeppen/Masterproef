function y = pert_kepler_sol(xspan,h)
y = [];
for i = xspan(1):h:xspan(2)
    sol = pert_keplerexact(i).';
    y = [y,sol];
end
end

function sol = pert_keplerexact(t)
global ep;
a=1+ep;
c=cos(a*t);
s=sin(a*t);
sol=[c,s,-a*s,a*c];
end
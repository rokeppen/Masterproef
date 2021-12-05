function fout = sympl3_simple(h,methode)
f=@simple;
J=@simple_jac;
sol=@simple_sol;

tspan = [0, 1000];
e=2;
y0 = [0,0,1,e];

if methode==1
[t,y] = sympl3_klassiek(f, tspan, y0 , J,h);
elseif methode==2
[t,y] = sympl3_Calvo(f, tspan, y0 , J,h);
elseif methode==3
[t,y] = sympl3_Marnix(f, tspan, y0 , J,h);
elseif methode==4
[t,y] = sympl3_Calvo_fixed(f, tspan, y0 , J,h);
end

for i=2:length(t)
  ysol=feval(sol,t(i));
  ff(i-1)=norm(y(i,:)-ysol,1);
end
fout=log10(max(ff));
end
function test_kepler(e)

 h=1/2;

f=@kepler;
J=@kepler_jac;
sol=@kepler_sol;

tspan = [0,1000];

y0 = [1-e,0,0,sqrt((1+e)/(1-e))];


[t,y] = sympl3_Calvo(f, tspan, y0 , J,h);
for i=2:length(t)
  ysol=feval(sol,t(i),e);
  fout(i-1)=(norm(y(i,:)-ysol,1));
end
fout;
max(fout)
log10(max(fout))
[t,y] = sympl3_Calvo_fixed(f, tspan, y0 , J,h);
for i=2:length(t)
  ysol=feval(sol,t(i),e);
  fout(i-1)=(norm(y(i,:)-ysol,1));
end
fout;
max(fout)
log10(max(fout))


end
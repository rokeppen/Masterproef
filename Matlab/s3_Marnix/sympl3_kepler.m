function fout=sympl3_kepler(h,methode);

f=@kepler;
J=@kepler_jac;
sol=@kepler_sol;

tspan = [0, 1000];
e=0.001;
y0 = [1-e,0,0,sqrt((1+e)/(1-e))];

if methode==1
[t,y] = sympl3_klassiek(f, tspan, y0 , J,h);
elseif methode==2
[t,y] = sympl3_Calvo(f, tspan, y0 , J,h);
elseif methode==3
[t,y] = sympl3_Marnix(f, tspan, y0 , J,h);
elseif methode==4
[t,y] = sympl3_Calvo_fixed(f, tspan, y0 , J,h);
end
% yend=y(end,:);
% ysol=feval(sol,tspan(2),e);
% fout=norm(yend-ysol,2);
for i=2:length(t)
  ysol=feval(sol,t(i),e);
  ff(i-1)=norm(y(i,:)-ysol,1);
end
fout=log10(max(ff))
end


   







        





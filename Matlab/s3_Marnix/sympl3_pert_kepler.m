function fout=sympl3_pert_kepler(h,methode);
global ep;
f=@pert_kepler;
J=@pert_kepler_jac;
sol=@pert_kepler_sol;

tspan = [0, 1000];
ep=0.001;
y0 = [1,0,0,1+ep];

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
% ysol=feval(sol,tspan(2));
% fout=norm(yend-ysol,2);
for i=2:length(t)
  ysol=feval(sol,t(i));
  ff(i-1)=norm(y(i,:)-ysol,2);
end
fout=log10(max(ff))
end



   







        





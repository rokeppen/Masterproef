function fout=sympl3_euler(h,methode);

global eulerbeta euleralpha;
%euleralpha=1+1/sqrt(1.51);
% eulerbeta=1-0.51/sqrt(1.51); 

f=@euler_prob;
J=@euler_jac;
sol=@euler_sol;

tspan = [0, 1000];
y0 = [0,1,1];

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
fout=max(ff)
end


   







        





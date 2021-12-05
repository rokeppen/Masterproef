function fout=sympl3_pendulum(h,methode);

f=@pendulum;
J=@pendulum_jac;
sol=@pendulum_sol;

tspan = [0, 100000];

y0 = [0,1.5];

if methode==1
[t,y] = sympl3_klassiek(f, tspan, y0 , J,h);
elseif methode==2
[t,y] = sympl3_Calvo(f, tspan, y0 , J,h);
elseif methode==3
[t,y] = sympl3_Marnix(f, tspan, y0 , J,h);
end
yend=y(end,:);
ysol=feval(sol,tspan(2));
fout=norm(yend(1)-ysol,2);

end


   







        





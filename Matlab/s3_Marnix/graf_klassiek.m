function grafiek
h=1;
err_klassiek=[];
err_Calvo=[];
hvals=[]
for i=1:4
  h=h/2;
  hvals=[hvals,log10(1000/h)];
  err_klassiek=[err_klassiek,log10(symplectisch3(h,1,@kepler,@kepler_jac,@kepler_sol))];
  err_Calvo=[err_Calvo,log10(symplectisch3(h,2,@kepler,@kepler_jac,@kepler_sol))];
end
plot(hvals,err_klassiek,'-c',hvals,err_Calvo,'-d')
xlabel('log10(steps)')
ylabel('log10(error)')
%title('van der Pol Equation, \mu = 1')
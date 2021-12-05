function grafiek_kepler

disp ('grafiek voor kepler probleem met s=3, klassiek en Calvo');
h=1;
err_klassiek=[];
err_Calvo=[];
hvals=[];
for i=1:4
  h=h/2;
  hvals=[hvals,log10(1000/h)];
  err_klassiek=[err_klassiek,log10(sympl3_kepler(h,1))];
  err_Calvo=[err_Calvo,log10(sympl3_kepler(h,2))];  % z=sqrt(-1)*h*(y(1)^2+y(2)^2)^(-3/2)
end
plot(hvals,err_klassiek,'-ob',hvals,err_Calvo,'-dr')
xlabel('log10(steps)')
ylabel('log10(error)')

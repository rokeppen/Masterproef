function grafiek_kepler
disp ('grafiek voor kepler probleem met s=3, Marnix met fac, en Calvo-fixed');
global fac;

h=1;
err_klassiek=[];
err_Calvo=[];
err_Marnix=[]; fac=1/2;
err_Calvo_fixed=[];
hvals=[];
for i=1:4
  h=h/2;
  hvals=[hvals,log10(1000/h)];
  err_klassiek=[err_klassiek,(sympl3_kepler(h,1))];
  err_Calvo=[err_Calvo,(sympl3_kepler(h,2))];  % z=sqrt(-1)*h*(y(1)^2+y(2)^2)^(-3/2)
  err_Marnix=[err_Marnix,(sympl3_kepler(h,3))];
  err_Calvo_fixed=[err_Calvo_fixed,(sympl3_kepler(h,4))];  % z=sqrt(-1)*h*(y(1)^2+y(2)^2)^(-3/2)
end
figure;
plot(hvals,err_klassiek,'-ob',hvals,err_Calvo,'-dr',hvals,err_Marnix,'-sm',hvals,err_Calvo_fixed,'-pg')
xlabel('log10(steps)')
ylabel('log10(error)')
title('Kepler problem');
legend('classic','Calvo','Van Daele','Calvo fixed')
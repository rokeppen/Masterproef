function grafiek_euler
disp ('grafiek voor euler probleem met s=3, klassiek, Calvo en Marnix en Calvo-fixed');
global fac
fac=1/2;
h=1;
err_klassiek=[];
err_Calvo=[];
err_Marnix=[];
err_Calvo_fixed=[];
hvals=[];
for i=1:4
  h=h/2;
  hvals=[hvals,log10(1000/h)];
  err_klassiek=[err_klassiek,log10(sympl3_euler(h,1))];
  err_Calvo=[err_Calvo,log10(sympl3_euler(h,2))];  % z=sqrt(-1)*h*(y(1)^2+y(2)^2)^(-3/2)
  err_Marnix=[err_Marnix,log10(sympl3_euler(h,3))];
  err_Calvo_fixed=[err_Calvo_fixed,log10(sympl3_euler(h,4))];
end
figure;
plot(hvals,err_klassiek,'-ob',hvals,err_Calvo,'-dr',hvals,err_Marnix,'-sm',hvals,err_Calvo_fixed,'-pg')
xlabel('log10(steps)')
ylabel('log10(error)')
title('Euler problem');
legend('classic','Calvo','Van Daele','Calvo fixed')
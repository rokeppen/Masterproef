function grafiek_euler_Marnix
disp ('grafiek voor Euler probleem met s=3, Marnix met fac');
global fac

h=1;
err_Marnix1=[];
err_Marnix2=[];
err_Marnix3=[];
err_Marnix4=[];
err_Marnix5=[];
err_Marnix6=[];
hvals=[];
for i=1:3
  h=h/2;
  hvals=[hvals,log10(1000/h)];
  fac=3
  err_Marnix1=[err_Marnix1,log10(sympl3_euler(h,3))];
  fac=0.5
  err_Marnix2=[err_Marnix2,log10(sympl3_euler(h,3))];
  fac=1.5
  err_Marnix3=[err_Marnix3,log10(sympl3_euler(h,3))];
  fac=2.5
  err_Marnix4=[err_Marnix4,log10(sympl3_euler(h,3))];
  fac=4
  err_Marnix5=[err_Marnix5,log10(sympl3_euler(h,3))];
  fac=0.25
  err_Marnix6=[err_Marnix6,log10(sympl3_euler(h,3))];

  
end
figure;
plot(hvals,err_Marnix1,'-ob',hvals,err_Marnix2,'-dr',hvals,err_Marnix3,'-sm',hvals,err_Marnix4,'--g',hvals,err_Marnix5,'-.k',...
     hvals,err_Marnix6,':pc')
xlabel('log10(steps)')
ylabel('log10(error)')
title('Euler problem');
legend('fac=3','fac=0.5','fac=1.5','fac=2.5','fac=4','fac=0.25');
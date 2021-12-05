function grafiek_simple_Marnix
global fac
h=1;
err_Marnix1=[];
hvals=[];
for i=1:3
  h=h/2;
  hvals=[hvals,log10(1000/h)];
  fac=2;
  err_Marnix1=[err_Marnix1,log10(sympl3_simple(h,3))]; 
end
figure;
plot(hvals,err_Marnix1,'-ob');
xlabel('log10(steps)')
ylabel('log10(error)')
title('Kepler problem');
legend('fac=3','fac=0.5','fac=1.5','fac=2.5','fac=4','fac=0.25');
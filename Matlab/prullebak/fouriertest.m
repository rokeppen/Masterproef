function fouriertest                    
T = 2*pi;
Fs = 1/T;
L = 100;
t = (0:L-1)*T;
X = kepler_exact(t);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
[psor,lsor] = findpeaks(P1,f,'SortStr','descend');
[psor(1),lsor(1)]
[psor(2),lsor(2)]
end

function sol = kepler_exact(t)
e=0.001;
E=t;
del=1;
while abs(del)>10^(-10)
  del=-(t-E+e.*sin(E))./(-1+e.*cos(E));
  E=E+del;
end
c=cos(E);
s=sin(E);
noem=1./(1-e*c);
sq=sqrt(1-e^2);
sol=[c-e,sq.*s,-s.*noem,sq.*noem.*c];
end
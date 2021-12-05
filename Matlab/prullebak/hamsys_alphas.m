function hamsys_alphas(inpute,xspan)
global e
e = inpute;
q0 = exact_sol([xspan(1),xspan(1)],1,'hamsys_exact').';
figure;
for i = 1:1
    h = 2^(-i);
    res = [];
    for alpha = 0:0.2:3
        q = RKs2(xspan,h,alpha,q0,@hamsystem,@(~,h) 1i*h);
        %q = RKs3(xspan,h,alpha,q0,@hamsystem,@(~,h) 1i*h);
        qe = exact_sol(xspan,h,'hamsys_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(0:0.2:3,res);
    hold on;
end
title('Euler: fout in functie van \alpha');
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

function y = hamsystem(q)
global e
k = 0.5;
a = e^2+k^2+1;
b = e^2-k^2-1;
y(1) = q(4);
y(2) = q(5);
%y(3) = -(a*q(2)^2+2*b*q(3)*q(2)+a*q(3)^2)*q(3)/4-q(3)*q(2)*(2*a*q(2)+2*b*q(3))/4+k^2/2*(q(2)-q(3))^3;
%y(4) = -(a*q(2)^2+2*b*q(3)*q(2)+a*q(3)^2)*q(2)/4-q(2)*q(3)*(2*a*q(3)+2*b*q(2))/4-k^2/2*(q(2)-q(3))^3;
y(3) = -(2*a*q(2)+2*b*q(3))/4+k^2/2*(q(2)-q(3))^3;
y(4) = -(2*a*q(3)+2*b*q(2))/4-k^2/2*(q(2)-q(3))^3;
end

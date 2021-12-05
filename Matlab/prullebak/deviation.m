function deviation(inpute,xspan)
global e
e = inpute;
q0 = exact_sol([xspan(1),xspan(1)],1,'simple_exact').';
alpha = corr_alpha(e);
range = alpha*0.9:alpha/100:alpha*1.1;
figure;
for i = 1:3
    h = 2^(-i);
    res = [];
    for a=range
        q = RKs2(xspan,h,a,q0,@simple_system,'testz2');
        qe = exact_sol(xspan,h,'simple_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
title("Simpel met \mu_1 = " + e + "*\mu_2: fout in functie van \alpha");
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

function dydt = simple_system(q)
global e
dydt(1) = q(4);
dydt(2) = q(5);
dydt(3) = -q(2);
dydt(4) = -q(3)*e^2;
end

function alpha = corr_alpha(k)
global oms
%oms = [-1,-6+3*k^2];
oms = [-k^2,3-6*k^2];
alpha = sqrt(oms(2)/oms(1));
end
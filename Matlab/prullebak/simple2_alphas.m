function simple2_alphas(inpute,xspan)
global e
e = inpute;
q0 = [1,1,1,e];
range = 1.2:0.2:3;
figure;
th = 0.288715183631118;
%th = NRtheta(-0.01,e,sqrt(3)/6);
for i = 1:1
    h = 0.1^i;
    res = [];
    for alpha=range
        q = RKs2(xspan,h,alpha,q0,@simple_system,'conszr',th,1i*h);
        qe = exact_sol(xspan,h,'simple2_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
title('Simpel met \mu_1 = 2*\mu_2: fout in functie van \alpha');
xlabel('\alpha');
ylabel('fout (1-norm)');
%legend('h=1/2','h=1/4','h=1/8');
end

function y = simple_system(q)
global e
y(1) = q(4);
y(2) = q(5);
y(3) = q(2);
y(4) = q(3)*e^2;
end
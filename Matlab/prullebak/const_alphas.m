function const_alphas(xspan)
q0 = [1];
%range = 0:0.2:3;
range = 1.2:0.2:3;
figure;
for i = 1:3
    h = 2^(-i);
    res = [];
    for alpha=range
        q = RKst(xspan,h,alpha,q0,@const_system,'consz');
        qe = ones(1,length(q));
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
title('Simpel met \mu_1 = 2*\mu_2: fout in functie van \alpha');
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

function dydt = const_system(q)
dydt(1) = q(2)-1;
end
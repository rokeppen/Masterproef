function simple_1d_alphas(inpute,xspan)
global e
e = inpute;
q0 = exact_sol([xspan(1),xspan(1)],1,'simple_1d_exact').';
range = 1.2:0.2:3;
figure;
for i = 1:2
    h = 2^(-i);
    res = [];
    for alpha=range
        q = RKs3(xspan,h,alpha,q0,@simple_system,@(~,h) h);
        qe = exact_sol(xspan,h,'simple_1d_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

function y = simple_system(q)
global e
y(1) = e*q(2);
%y(1) = e*exp(e*q(1));
end
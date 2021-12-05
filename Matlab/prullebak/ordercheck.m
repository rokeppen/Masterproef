function rico = ordercheck
e = 0.001;
global ep;
ep = 0.001;
xspan = [0,100];
%q0 = [1-e,0,0,sqrt((1+e)/(1-e))];
q0 = [1,0,0,1+ep];
alpha = 1/2;
x = [];
res = [];
for i = 1:4
    h = 2^(-i);
    %q = RKs2_kepler(xspan,h,alpha,q0,@keplersystem);
    %qe = kepler_sol(xspan,h,e);
    q = RKs2_pert_kepler(xspan,h,alpha,q0,@pert_keplersystem);
    qe = pert_kepler_sol(xspan,h);
    res = [res,log10(max(norm(q-qe,1)))];
    x = [x,log10(h)];
end
figure;
plot(x,res);
title('Log-log-diagram');
xlabel('log(h)');
ylabel('log(fout)');
rico = (res(end)-res(1))/(x(end)-x(1));
end

function y = keplersystem(q1,q2,p1,p2)
y(1) = p1;
y(2) = p2;
y(3) = -q1/(q1^2+q2^2)^(3/2);
y(4) = -q2/(q1^2+q2^2)^(3/2);
end

function dydt = pert_keplersystem(q1,q2,p1,p2)
global ep
noem=sqrt(q1^2+q2^2);
noem3=noem^3;
noem5=noem^5;
dydt = [p1;
        p2;
        -q1/noem3-(2*ep+ep^2)*q1/noem5;
        -q2/noem3-(2*ep+ep^2)*q2/noem5];
end
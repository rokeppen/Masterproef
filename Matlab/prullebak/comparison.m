function comparison
xspan = [0,10];
global e
e = 3;
q0 = exact_sol([xspan(1),xspan(1)],1,'simple_exact').';
range = [];
res2 = [];
res3 = [];
rest = [];
for i = 1:4
    h = 2^(-i);
    range = [range,log10(h)];
    q2 = RKs2(xspan,h,corr_alpha(e),q0,@simple_system,'testz');
    q3 = RKs3(xspan,h,0,q0,@simple_system,'consz');
    %qt = RKs2_calvo(xspan,h,q0,@simple_system);
    qe = exact_sol(xspan,h,'simple_exact');
    res2 = [res2,log10(max(norm(q2-qe,1)))];
    res3 = [res3,log10(max(norm(q3-qe,1)))];
    %rest = [rest,log10(max(norm(q-qe,1)))];
end
figure;
plot(range,res2);
hold on;
plot(range,res3);
%hold on;
%plot(range,rest);
hold off;
title('Fout in functie van h');
xlabel('log10(h)');
ylabel('log10(fout)');
legend('RKs2','RKs3','Conventional');
end

function dydt = simple_system(q)
%q = [t,y1,y2,y3,y4]
global e
dydt(1) = q(4); %cos(q(1));
dydt(2) = q(5); %e*cos(e*q(1));
dydt(3) = -q(2); %-sin(q(1));
dydt(4) = -q(3)*e^2; %-sin(e*q(1))*e^2;
end

function alpha = corr_alpha(k)
global oms
oms = [-k^2,3-6*k^2];
alpha = sqrt(oms(2)/oms(1));
end
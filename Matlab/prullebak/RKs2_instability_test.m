function RKs2_instability_test(max)
q0 = 1;
figure;
global h
h = 0.1;
theta = sqrt(3)/6;
q = RKs2([0,max],h,1,q0,@simple_system,@(~,~) pi*1i/(4*theta) - 1);
qe = exact_sol([0,max],h,'simple_1d_test');
plot(0:h:max,q);
legend('RKs2');
figure;
plot(0:h:max,qe);
legend('exact');
%legend('RKs2','Exacte oplossing');
end

function y = simple_system(q)
global h
global omega
omega = 2+6*1i; % instabiel, maar numeriek stabiel
%omega = -3-10*1i; % stabiel, maar numeriek instabiel
y(1) = omega/h*q(2);
end
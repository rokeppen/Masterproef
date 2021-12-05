function pendulum_alphas(inita)
global a
a = inita;
xspan = [0,100000];
q0 = [0,1.5];
range=0:0.2:3;
res = [];
figure;
h=1/2;
for alpha=range
    q = RKs2_pendulum(xspan,h,alpha,q0,@pendulumsystem);
    qe = -0.595399559;
    res = [res,log10(q(1,end)-qe)];
end
title('Pendulum: fout in functie van \alpha');
xlabel('\alpha');
ylabel('fout (1-norm)');
end

function dydt = pendulumsystem(q)
global a
dydt = [q(2);
        -a*sin(q(1))];
end
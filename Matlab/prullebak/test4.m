function test4(xspan)
global e
e = 2;
q0 = [0,0,1,e];
%figure;
testje = 10000;
th = 0.288715183631118; %opt
0.288721630783035;
for ths = th*0.9999999:th/100000000:th*1.0000001
    for i = 1:1
        h = 0.1^i;
        res = [];
        for alpha=2%range
            z=1i*h;
            q = RKs2(xspan,h,alpha,q0,@simple_system,'consz',ths,z);
            qe = exact_sol(xspan,h,'simple_exact');
            res = [res,log10(max(norm(q-qe,1)))];
        end
        %plot(range,res);
        %hold on;
    end
    t = max(norm(q-qe,1));
    if testje > t
        testje = t;
        thopt = ths;
    end
end
thopt
%title('Simpel met \mu_1 = 2*\mu_2: fout in functie van \alpha');
%xlabel('\alpha');
%ylabel('fout (1-norm)');
%legend('h=1/2','h=1/4','h=1/8');
end

function y = simple_system(q)
global e
y(1) = q(4);
y(2) = q(5);
y(3) = -q(2);
y(4) = -q(3)*e^2;
end
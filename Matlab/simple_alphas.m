% plot voor een range van verhoudingen alpha de fout van de methode,
% toegepast op het lineair autonoom probleem met k_1=i en k_2=e
%   @param inpute: de waarde voor k_2
%   @param range: het interval waarin alpha wordt bekeken
%   @param xspan: het interval waarover geïntegreerd wordt
%   @param scheme: de te gebruiken numerieke methode
function simple_alphas(inpute,range,xspan,scheme)
global e
e = inpute;
q0 = exact_sol([xspan(1),xspan(1)],1,'simple_exact').';
figure;
for i = 1:3
    h = 2^(-i);
    res = [];
    for alpha = range
        q = scheme(xspan,h,sqrt(alpha),q0,@simple_system,@(~,h) e*1i*h);
        qe = exact_sol(xspan,h,'simple_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

% stelsel voor het lineair autonoom probleem 
function dydt = simple_system(q)
%q = [t,y1,y2,y3,y4]
global e
dydt(1) = q(4); %cos(q(1));
dydt(2) = q(5); %e*cos(e*q(1));
dydt(3) = -q(2); %-sin(q(1));
dydt(4) = -q(3)*e^2; %-sin(e*q(1))*e^2;
end
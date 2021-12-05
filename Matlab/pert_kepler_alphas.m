% plot voor een range van verhoudingen alpha de fout van de methode,
% toegepast op het geperturbeerde tweelichamenprobleem van kepler
%   @param range: het interval waarin alpha wordt bekeken
%   @param xspan: het interval waarover geïntegreerd wordt
%   @param scheme: de te gebruiken numerieke methode
function pert_kepler_alphas(range,xspan,scheme)
global e
e = 0.001; % perturbatiefactor
q0 = exact_sol([xspan(1),xspan(1)],1,'pert_kepler_exact').';
figure;
for i = 1:3
    h = 2^(-i);
    res = [];
    for alpha = range
        q = scheme(xspan,h,sqrt(alpha),q0,@pert_keplersystem,@(~,h) (1+e)*1i*h);
        qe = exact_sol(xspan,h,'pert_kepler_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

% stelsel voor het geperturbeerde keplerprobleem
function dydt = pert_keplersystem(q)
global e
n = sqrt(q(2)^2+q(3)^2);
dydt = [q(4),q(5),-q(2)/n^3-(2*e+e^2)*q(2)/n^5,-q(3)/n^3-(2*e+e^2)*q(3)/n^5];
end
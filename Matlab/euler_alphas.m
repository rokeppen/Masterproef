% plot voor een range van verhoudingen alpha de fout van de methode,
% toegepast op het probleem van euler
%   @param range: het interval waarin alpha wordt bekeken
%   @param xspan: het interval waarover geïntegreerd wordt
%   @param scheme: de te gebruiken numerieke methode
function euler_alphas(range,xspan,scheme)
q0 = exact_sol([xspan(1),xspan(1)],1,'euler_exact').';
figure;
for i = 1:3
    h = 2^(-i);
    res = [];
    for alpha = range
        q = scheme(xspan,h,sqrt(alpha),q0,@eulersystem,@(~,h) 1i*h*2*pi/7.45056320933095);
        qe = exact_sol(xspan,h,'euler_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

% stelsel voor het probleem van euler
function y = eulersystem(q)
a = 1+1/sqrt(1.51);
b = 1-0.51/sqrt(1.51);
y(1) = (a-b)*q(3)*q(4);
y(2) = (1-a)*q(2)*q(4);
y(3) = (b-1)*q(2)*q(3);
end
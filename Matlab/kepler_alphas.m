% plot voor een range van verhoudingen alpha de fout van de methode,
% toegepast op het tweelichamenprobleem van kepler
%   @param range: het interval waarin alpha wordt bekeken
%   @param xspan: het interval waarover geïntegreerd wordt
%   @param scheme: de te gebruiken numerieke methode
function kepler_alphas(range,xspan,scheme)
global e
e = 0.001; % excentriciteit
q0 = exact_sol([xspan(1),xspan(1)],1,'kepler_exact').';
figure;
for i = 1:3
    h = 2^(-i);
    res = [];
    for alpha = range
        % variabele tweede frequentie:
        %q = scheme(xspan,h,sqrt(alpha),q0,@kepler_system,@(y,h) 1i*h*(y(1)^2+y(2)^2)^(-3/2));
        q = scheme(xspan,h,sqrt(alpha),q0,@kepler_system,@(~,h) 1i*h);
        qe = exact_sol(xspan,h,'kepler_exact');
        res = [res,log10(max(norm(q-qe,1)))];
    end
    plot(range,res);
    hold on;
end
xlabel('\alpha');
ylabel('fout (1-norm)');
legend('h=1/2','h=1/4','h=1/8');
end

% stelsel voor het keplerprobleem
function y = kepler_system(q)
y(1) = q(4);
y(2) = q(5);
y(3) = -q(2)/(q(2)^2+q(3)^2)^(3/2);
y(4) = -q(3)/(q(2)^2+q(3)^2)^(3/2);
end
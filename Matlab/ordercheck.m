% log-log-diagram ter illustratie van de mogelijke ordeverhoging bij een
% gepaste keuze van de fitting space
%   @param scheme: de te gebruiken numerieke methode
%   @param inpute: de tweede te gebruiken frequentie voor het lineair
%                  autonoom probleem
%   @param range: de range voor i ter bepaling van de stapgroottes h=2^(-i)
function ordercheck(scheme,inpute,range)
global e;
global oms;
e = inpute;
q0 = exact_sol([0,0],1,'simple_exact').';
figure;
alphas = [corr_alpha(e,func2str(scheme)), e];
for j = 1:2
    alpha = alphas(j);
    res = [];
    x = [];
    for i = range
        h = 2^(-i);
        xspan = [0,h];
        q = scheme(xspan,h,alpha,q0,@simple_system,@(~,h) sqrt(oms(1))*h);
        qe = exact_sol(xspan,h,'simple_exact');
        res = [res,log10(max(norm(q-qe,1)))];
        x = [x,log10(h)];
    end
    plot(x,res);
    hold on;
    % fout = Ch^(p+1)+O(h^(p+2))
    % => log(fout) ~ log(C)+(p+1)*log(h)
    % => p = rico-1
    order = (res(end)-res(1))/(x(end)-x(1))-1
end
xlabel('log(h)');
ylabel('log(fout)');
legend({'optimale frequenties','niet-optimale frequenties'},'Location','northwest');
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

% bepaal de optimale verhouding alpha afgeleid uit de foutterm
function alpha = corr_alpha(k,scheme)
global oms
if scheme == "RKs2"
    oms = [-1,-6+3*k^2];
    %oms = [-k^2,3-6*k^2]; % optie 2
elseif scheme == "RKs3"
    oms = [-1,10*(k^2-k^4)/(2-3*k^2)];
    %oms = [-k^2,10*(k^2-1)/(2*k^2-3)]; % optie 2
else
    error("Niet geïmplementeerd schema");
end
alpha = sqrt(oms(2)/oms(1));
end
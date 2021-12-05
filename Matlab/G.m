% hulpfunctie bij de bepaling van theta
function y = G(a,b,t2)
if abs(a-b) < 0.01 
    y = (eta(-1,a/4)-eta(0,a/4))/a/t2/eta(0,t2*a);
else   
    y = (eta(0,a/4)-eta(0,b/4))/(eta(-1,a*t2)-eta(-1,b*t2));
end
end
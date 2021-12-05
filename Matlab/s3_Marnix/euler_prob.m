function dydt = euler_prob(y)
global euleralpha eulerbeta;
% euleralpha=1+1/sqrt(1.51);
% eulerbeta=1-0.51/sqrt(1.51);
dydt = [(euleralpha-eulerbeta)*y(2)*y(3);
        (1-euleralpha)*y(1)*y(3);
        (eulerbeta-1)*y(1)*y(2)];
end
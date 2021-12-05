function dydt = euler(y)
global euleralpha eulerbeta;
dydt = [(euleralpha-eulerbeta)*y(2)*y(3);
        (1-euleralpha)*y(1)*y(3);
        (eulerbeta-1)*y(1)*y(2)];
end
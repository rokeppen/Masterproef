function dydt = pendulum(yp)
global a_pend
dydt = [y(2);
        -a_pend*sin(y(1))];
end
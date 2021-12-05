function dfdy = pendulum_jac(y)
global a_pend
a_pend=5;
dfdy= [0,1;
       -a_pend*cos(y(1)),0];
end


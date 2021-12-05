function dfdy = euler_jac(y)
global eulerbeta euleralpha;
euleralpha=1+1/sqrt(1.51);
eulerbeta=1-0.51/sqrt(1.51);    
dfdy= [0,(euleralpha-eulerbeta)*y(3),(euleralpha-eulerbeta)*y(2);
       (1-euleralpha)*y(3),0,(1-euleralpha)*y(1);
       (eulerbeta-1)*y(2),(eulerbeta-1)*y(1),0];
end
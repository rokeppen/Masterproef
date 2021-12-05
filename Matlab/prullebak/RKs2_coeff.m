function RKs2_coeff
syms gam th lam b1 h mu1 z;
h = z/mu1;
c = [1/2-th, 1/2+th];
b = [b1, b1];
a = [gam*b1/2, gam*b1/2+lam; gam*b1/2-lam, gam*b1/2];
fun = [Lexp(a(1,:), c, c(1), gam, h, mu1), Lexp(a(2,:), c, c(2), gam, h, mu1), Lexp(b, c, 1, 1, h, mu1)];
sol = solve(fun,[gam, lam, b1]);
gamsol = simplify(sol.gam)
lamsol = simplify(sol.lam)
b1sol = simplify(sol.b1)
%thsol = simplify(sol.th)
end

function eq = Lexp(a, c, ci, gam, h, mu)
syms x;
eq = exp(mu*(x+h*ci))-gam*exp(mu*x)-h*a(1)*mu*exp(mu*(x+c(1)*h))-h*a(2)*mu*exp(mu*(x+c(2)*h)) == 0;
end
function RKs3_coeff(opt)
syms g1 g2 th b1 b2 a2 a3 a4 h m z;
h = z/m;
c = [1/2-th, 1/2, 1/2+th];
g = [g1, g2, g1];
b = [b1, b2, b1];
a = [g1*b1/2, g1*b2/2-a2, g1*b1/2-a3; g2*b1/2-a4, g2*b2/2, g2*b1/2+a4; g1*b1/2+a3, g1*b2/2+a2, g1*b1/2];
% extra voorwaarde voor symmetrie/symplecticiteit
eq1 = b1*a2/g1+b2*a4/g2 == 0;
% eerste frequentie
eq2 = L(a(1,:), c, c(1), g(1), h, @f1);
eq3 = L(a(2,:), c, c(2), g(2), h, @f1);
eq4 = L(a(3,:), c, c(3), g(3), h, @f1);
eq5 = L(b, c, 1, 1, h, @f1);
if opt == 1
    % Optie 1: f(x) = x voor externe trap
    eqs = [L(b, c, 1, 1, h, @(x) x), L(b, c, 1, 1, h, @(x) x^2)];
elseif opt == 2
    %Optie 2: RKs3 (constante gamma's)
    eqs = [L(a(1,:), c, c(1), g(1), h, @(x) x^0), L(a(2,:), c, c(2), g(2), h, @(x) x^0)];
elseif opt == 3
    %Optie 3: f(x) = x*exp(x) voor externe trap
    eqs = L(b, c, 1, 1, h, @(x) x*exp(m*x));
elseif opt == 4
    %Optie 4: Extra frequentie op de externe trap
    syms n;
    eqs = L(b, c, 1, 1, h, @(x) exp(n*x));
else
    eqs = [];
end
sol = solve([eq1, eq2, eq3, eq4, eq5, eqs],[g1, g2, b1, b2, a2, a3, a4]);
g1sol = simplify(sol.g1)
g2sol = simplify(sol.g2)
b1sol = sol.b1 %simplify(sol.b1)
b2sol = simplify(sol.b2)
a2sol = simplify(sol.a2)
a3sol = simplify(sol.a3)
a4sol = simplify(sol.a4)
end

function eq = L(a, c, ci, gam, h, f)
syms x;
eq = f(x+h*ci)-gam*f(x)-h*a(1)*diff(f(x+c(1)*h))-h*a(2)*diff(f(x+c(2)*h))-h*a(3)*diff(f(x+c(3)*h)) == 0;
end

function y = f1(x)
syms m;
y = exp(m*x);
end
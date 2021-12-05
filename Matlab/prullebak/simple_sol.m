function y = simple_sol(xspan,h)
y = [];
for i = xspan(1):h:xspan(2)
    sol = simple_exact(i).';
    y = [y,sol];
end
end

function sol = simple_exact(t)
global a;
sol=[sin(t),sin(a*t),cos(t),a*cos(a*t)];
end
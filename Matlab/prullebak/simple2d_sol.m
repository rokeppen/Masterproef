function y = simple2d_sol(xspan,h)
y = [];
for i = xspan(1):h:xspan(2)
    sol = simple_exact(i).';
    y = [y,sol];
end
end


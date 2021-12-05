function y = euler_sol(xspan,h)
y = [];
for i = xspan(1):h:xspan(2)
    sol = eulerexact(i).';
    y = [y,sol];
end
end


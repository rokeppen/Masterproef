% functies van Ixaru
%   @param index: index van de functie, tussen -1 en oneindig
%   @param z: argument van de functie
function y = eta(index, z)
if index == -1
    if z <= -0.1
        y = cos(sqrt(abs(z)));
    elseif abs(z) < 0.1
        y = z^15/265252859812191058636308480000000 + z^14/304888344611713860501504000000 + z^13/403291461126605635584000000 + z^12/620448401733239439360000 + z^11/1124000727777607680000 + z^10/2432902008176640000 + z^9/6402373705728000 + z^8/20922789888000 + z^7/87178291200 + z^6/479001600 + z^5/3628800 + z^4/40320 + z^3/720 + z^2/24 + z/2 + 1;
    else
        y = cosh(sqrt(z));
    end
elseif index == 0
    if z <= -0.1
        y = sin(sqrt(abs(z)))/sqrt(abs(z));
    elseif abs(z) < 0.1
        y = z^15/8222838654177922817725562880000000 + z^14/8841761993739701954543616000000 + z^13/10888869450418352160768000000 + z^12/15511210043330985984000000 + z^11/25852016738884976640000 + z^10/51090942171709440000 + z^9/121645100408832000 + z^8/355687428096000 + z^7/1307674368000 + z^6/6227020800 + z^5/39916800 + z^4/362880 + z^3/5040 + z^2/120 + z/6 + 1;
    else
        y = sinh(sqrt(z))/sqrt(z);
    end
elseif z == 0
    y = 1/factorial(factorial(2*index+1));
else
    y = (eta(index - 2, z) - (2*index - 1)*eta(index - 1, z))/z;
end
end
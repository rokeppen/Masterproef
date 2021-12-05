function dydt = pert_kepler(y)
global ep
noem=(y(1)^2+y(2)^2)^(1/2);
noem3=noem^3;
noem5=noem^5;
dydt = [y(3);
        y(4);
        -y(1)/noem3-(2*ep+ep^2)*y(1)/noem5;
        -y(2)/noem3-(2*ep+ep^2)*y(2)/noem5];
end
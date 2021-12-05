function dydt = twodim_oscil(y)
noem=(y(1)^2+y(2)^2)^(1/2);
noem3=noem^3;
dydt = [y(3);
        y(4);
        -y(1)/noem3;
        -y(2)/noem3];
end
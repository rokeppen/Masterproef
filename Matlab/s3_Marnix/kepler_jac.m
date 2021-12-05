function dfdy = kepler_jac(y)
noem=(y(1)^2+y(2)^2)^(1/2);
noem3=noem^3;
noem5=noem^5;
dfdy= [0,0,1,0;0,0,0,1;
       -1/noem3+3*y(1)^2/noem5,3*y(1)*y(2)/noem5,0,0;
       3*y(1)*y(2)/noem5,-1/noem3+3*y(2)^2/noem5,0,0];
end


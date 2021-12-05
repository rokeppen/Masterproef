function dfdy = pert_kepler_jac(y)
global ep;
a=(2*ep+ep^2)/3;
noem=(y(1)^2+y(2)^2)^(1/2);
noem3=noem^3;
noem5=noem^5;
noem7=noem^7;
dfdy= [0,0,1,0;0,0,0,1;
       -1/noem3+3*y(1)^2/noem5-a*y(1)/noem5+5*a*y(1)^2/noem7,3*y(1)*y(2)/noem5+5*a*y(1)*y(2)/noem7,0,0;
       3*y(1)*y(2)/noem5+5*a*y(1)*y(2)/noem7,-1/noem3+3*y(2)^2/noem5-a*y(2)/noem5+5*a*y(2)^2/noem7,0,0];
end


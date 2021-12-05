function y = test
th = NRtheta4_bis(z,z2,[sqrt((15+2*sqrt(30))/35)/2, sqrt((15-2*sqrt(30))/35)/2]);
a11 = @(z,z2,theta1,theta2) (-z*sinh(z*theta2)*(cosh(z2*theta1)*cosh(z*theta2) - cosh(z2*theta2)*cosh(z*theta1))*cosh(z2/2) + z*cosh(z*theta2)*(sinh(z2*theta2)*sinh(z*theta1) - sinh(z*theta2)*sinh(z2*theta1))*sinh(z2/2) + sinh(z2*theta2)*z2*(cosh(z2*theta1)*cosh(z*theta2) - cosh(z2*theta2)*cosh(z*theta1))*cosh(z/2) - cosh(z2*theta2)*z2*(sinh(z2*theta2)*sinh(z*theta1) - sinh(z*theta2)*sinh(z2*theta1))*sinh(z/2) + ((2*z2*cosh(z*theta1)^2 - z2)*sinh(z2*theta2) - sinh(z*theta2)*(z*cosh(z2*theta1)*cosh(z*theta1) + sinh(z*theta1)*sinh(z2*theta1)*z2))*cosh(z2*theta2) - ((z*sinh(z*theta1)*sinh(z2*theta1) + z2*cosh(z*theta1)*cosh(z2*theta1))*sinh(z2*theta2) - 2*z*(cosh(z2*theta1)^2 - 1/2)*sinh(z*theta2))*cosh(z*theta2))/2*z2*z*(cosh(z2*theta1)*cosh(z*theta2) - cosh(z2*theta2)*cosh(z*theta1))*(sinh(z2*theta2)*sinh(z*theta1) - sinh(z*theta2)*sinh(z2*theta1));
a11s = @(z,z2,theta1,theta2) ;
end


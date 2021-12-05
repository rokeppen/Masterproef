% bepaling van de theta's voor RKs4
%   @param z1: waarde voor mu1*h
%   @param z2: waarde voor mu2*h
function y = NRtheta4(z1,z2)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
% gaussische integratiepunten bij het klassieke geval
prev = [sqrt((15+2*sqrt(30))/35)/2, sqrt((15-2*sqrt(30))/35)/2];
% threshold vanaf wanneer reeksontwikkelingen dienen gebruikt te worden
thresh = 1.2;
if z1 == 0 && z2 == 0
    y = prev;
elseif abs(z1) < thresh && abs(z2) < thresh
    z = z1;
    y(1) = sqrt(525 + 70*sqrt(30))/70 + ((18289*sqrt(30) + 305094)*(z^2 + z2^2))/(14112*sqrt(525 + 70*sqrt(30))*(82*sqrt(30) + 16813)) - (390784237170926385557973829*sqrt(30) + 2255680480321825145698558710)*sqrt(525 + 70*sqrt(30))*(z^4 + z2^4)/(157228731437834403203936225642496000*sqrt(30) + 870339921440433446260756785184128000) - 4.611063613*10^(-8)*z2^2*z^2 + (264791285795758401521259928894835496973586021377162084626316485215287641771108149709907*sqrt(30) + 1450318634900293761183179222213988800243229653278786115275125509271640207716025863815280)*sqrt(525 + 70*sqrt(30))*(z^6 + z2^6)/(79708094356773760148542339222050212306336998333939553262863315961352549787772001660134923089920000*sqrt(30) + 436579093531717448423931212999765249176438101886253199133678699129500237236379566711945034332160000) + 1.777228373*10^(-10)*z2^2*z^4 + 1.777228373*10^(-10)*z2^4*z^2;
    y(2) = sqrt(525 - 70*sqrt(30))/70 + (18289*sqrt(30) + 305094)*(z^2 + z2^2)/((5 + 3*sqrt(30))^(3/2)*sqrt(-3 + sqrt(30))*(165312*sqrt(30) + 33895008)) - sqrt(5 + 3*sqrt(30))*sqrt(-3 + sqrt(30))*(991*sqrt(30) + 2448)*(z^4 + z2^4)/383360947200 + 2.340907267*10^(-7)*z2^2*z^2 + sqrt(5 + 3*sqrt(30))*sqrt(-3 + sqrt(30))*(5702597*sqrt(30) + 18716356)*(z^6 + z2^6)/1230772653766656000 - 3.872086760*10^(-10)*z2^2*z^4 - 3.872086760*10^(-10)*z2^4*z^2;
elseif abs(z1) < 0.4 && z1 ~= 0
    z = z1;
    fun = @(t) [-(-cosh(z2*t(1))*t(2)^2*z2/2 + t(1)^2*cosh(z2*t(2))*z2/2 + t(2)^2*sinh(z2/2) - t(1)^2*sinh(z2/2) - cosh(z2*t(2))*z2/24 + cosh(z2*t(1))*z2/24)*z^2/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2) - (-cosh(z2*t(1))*t(2)^4*z2/24 + t(1)^4*cosh(z2*t(2))*z2/24 + t(2)^4*sinh(z2/2)/12 - t(1)^4*sinh(z2/2)/12 - cosh(z2*t(2))*z2/1920 + cosh(z2*t(1))*z2/1920 + (-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^2*cosh(z2*t(2))/2 - cosh(z2*t(1))*t(2)^2/2)/(24*(cosh(z2*t(2)) - cosh(z2*t(1)))))*z^4/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2),
        -(-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2) + (-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^2*cosh(z2*t(2))/2 - cosh(z2*t(1))*t(2)^2/2)*z^2/((cosh(z2*t(2)) - cosh(z2*t(1)))^2*z2) - (t(2)^4*sinh(z2/2)*t(1)^2 - t(1)^4*sinh(z2/2)*t(2)^2 - t(1)^2*cosh(z2*t(2))*z2/160 + cosh(z2*t(1))*t(2)^2*z2/160 - cosh(z2*t(1))*t(2)^4*z2/24 + t(1)^4*cosh(z2*t(2))*z2/24 - (-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^4*cosh(z2*t(2))/24 - cosh(z2*t(1))*t(2)^4/24)/(cosh(z2*t(2)) - cosh(z2*t(1))) + ((-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^2*cosh(z2*t(2)) - cosh(z2*t(1))*t(2)^2)*(t(1)^2*cosh(z2*t(2))/2 - cosh(z2*t(1))*t(2)^2/2))/(2*(cosh(z2*t(2)) - cosh(z2*t(1)))^2))*z^4/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2)];
    y = fsolve(fun,prev,opt);
elseif abs(z2) < 0.4 && z2 ~= 0
    z = z2;
    z2 = z1;
    fun = @(t) [-(-cosh(z2*t(1))*t(2)^2*z2/2 + t(1)^2*cosh(z2*t(2))*z2/2 + t(2)^2*sinh(z2/2) - t(1)^2*sinh(z2/2) - cosh(z2*t(2))*z2/24 + cosh(z2*t(1))*z2/24)*z^2/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2) - (-cosh(z2*t(1))*t(2)^4*z2/24 + t(1)^4*cosh(z2*t(2))*z2/24 + t(2)^4*sinh(z2/2)/12 - t(1)^4*sinh(z2/2)/12 - cosh(z2*t(2))*z2/1920 + cosh(z2*t(1))*z2/1920 + (-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^2*cosh(z2*t(2))/2 - cosh(z2*t(1))*t(2)^2/2)/(24*(cosh(z2*t(2)) - cosh(z2*t(1)))))*z^4/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2),
        -(-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2) + (-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^2*cosh(z2*t(2))/2 - cosh(z2*t(1))*t(2)^2/2)*z^2/((cosh(z2*t(2)) - cosh(z2*t(1)))^2*z2) - (t(2)^4*sinh(z2/2)*t(1)^2 - t(1)^4*sinh(z2/2)*t(2)^2 - t(1)^2*cosh(z2*t(2))*z2/160 + cosh(z2*t(1))*t(2)^2*z2/160 - cosh(z2*t(1))*t(2)^4*z2/24 + t(1)^4*cosh(z2*t(2))*z2/24 - (-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^4*cosh(z2*t(2))/24 - cosh(z2*t(1))*t(2)^4/24)/(cosh(z2*t(2)) - cosh(z2*t(1))) + ((-12*t(1)^2*cosh(z2*t(2))*z2 + 12*cosh(z2*t(1))*t(2)^2*z2 + 24*t(1)^2*sinh(z2/2) - 24*t(2)^2*sinh(z2/2) + cosh(z2*t(2))*z2 - cosh(z2*t(1))*z2)*(t(1)^2*cosh(z2*t(2)) - cosh(z2*t(1))*t(2)^2)*(t(1)^2*cosh(z2*t(2))/2 - cosh(z2*t(1))*t(2)^2/2))/(2*(cosh(z2*t(2)) - cosh(z2*t(1)))^2))*z^4/((cosh(z2*t(2)) - cosh(z2*t(1)))*z2)];
    y = fsolve(fun,prev,opt);
elseif z1 == 0 || z2 == 0
    if z1 == 0
        z = z2;
    else
        z = z1;
    end
    fun = @(t) [240*t(1)^2*t(2)^2-20*t(1)^2-20*t(2)^2+3;
        -12*z*((-t(2)^2+1/12)*cosh(z*t(1))+cosh(z*t(2))*(t(1)^2-1/12))*cosh(z/2)+(12*t(1)^2-12*t(2)^2)*sinh(z)];
    y = fsolve(fun,prev,opt);
elseif abs(z2) == abs(z1)
    fun = @(t) [2*b1bis(z1,t(1),t(2))+2*b2bis(z1,t(1),t(2))-1; 
        24*b1bis(z1,t(1),t(2))*t(1)^2+24*b2bis(z1,t(1),t(2))*t(2)^2-1];
    y = fsolve(fun,prev,opt);
else
    fun = @(t) [2*b1(z1^2,z2^2,t(1)^2,t(2)^2)+2*b2(z1^2,z2^2,t(1)^2,t(2)^2)-1; 
        24*b1(z1^2,z2^2,t(1)^2,t(2)^2)*t(1)^2+24*b2(z1^2,z2^2,t(1)^2,t(2)^2)*t(2)^2-1];
    y = fsolve(fun,prev,opt);
end
end

% b1-coëfficiënt bij RKs4
function y = b1(z,z2,t1,t2)
y = (eta(-1,z2*t2)*eta(0,z/4)-eta(-1,z*t2)*eta(0,z2/4))/2/(eta(-1,z2*t2)*eta(-1,z*t1)-eta(-1,z2*t1)*eta(-1,z*t2));
end

% b2-coëfficiënt bij RKs4
function y = b2(z,z2,t1,t2)
y = b1(z,z2,t2,t1);
end

% afgeleide van b1-coëfficiënt bij RKs4
function y = b1bis(z,theta1,theta2)
y = ((-2*theta2*(z*cosh(z*theta1)*theta1 - sinh(z*theta1))*cosh(z*theta2)^2 + 2*sinh(z*theta2)*(sinh(z*theta1)*theta2^2*z - theta1*cosh(z*theta1))*cosh(z*theta2) + 2*cosh(z*theta1)*theta1*theta2*z)*sinh(z/2) + z*cosh(z/2)*cosh(z*theta2)*(theta1*cosh(z*theta1)*sinh(z*theta2) - theta2*sinh(z*theta1)*cosh(z*theta2)))/(2*(-2*theta1*(cosh(z*theta1)^2 - 1/2)*theta2*cosh(z*theta2)^2 + cosh(z*theta1)*sinh(z*theta1)*sinh(z*theta2)*(theta1^2 + theta2^2)*cosh(z*theta2) + cosh(z*theta1)^2*theta1*theta2)*z^2);
end

% afgeleide van b2-coëfficiënt bij RKs4
function y = b2bis(z,t1,t2)
y = b1bis(z,t2,t1);
end
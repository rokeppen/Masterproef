% bepaling van de theta's voor RKs3
%   @param z: waarde voor mu1*h
%   @param z2: waarde voor mu2*h
function th = NRtheta3(z,z2)
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
% threshold vanaf wanneer reeksontwikkelingen dienen gebruikt te worden
d = 0.01;
if (z2==0) && (abs(z)>d)
    th = 1/sqrt(z)*acosh((1+eta(0,z)-2*eta(0,z/4))/2/(eta(0,z/4)-1));
elseif (z==0) && (abs(z2)>d)
    th=sqrt(15)/10;
    del=1;
    iter=0;
    while (abs(del)>10^(-10)) && (iter<10)
       th2=th^2;
       del= -((eta(-1,z2*th2)-1)/12/th2+1-eta(0,z2/4))/(z2/12*eta(0,z2*th2)-(eta(-1,z2*th2)-1)/6/th2)*th;
       th=th+del;
       iter=iter+1;
    end
elseif (abs(z)>d) && (abs(z2)>d)
    if abs(z2-4*z)<0.0001
       z1=z2;
    else
       z1=z;
    end
    th = fsolve(@(t) G(z1,z2,t^2)-G(z,4*z,t^2),sqrt(15)/10,opt);
else
    th = sqrt(15)/10*(1+(1/2100+(-131/105840000+(13487/4889808000000+ ...
           (-1175117/320380220160000000-505147/91537205760000000000*z2)*z2)*z2)*z2)*z2 ...
            +(1/420+(-17/21168000+(-1153/1955923200000+(-21517/16019011008000000+137477/4660075929600000000*z2)*z2)*z2)*z2...
            +(-17/784000+(-330733/1955923200000+(223250821/320380220160000000-43658273/25630417612800000000*z2)*z2)*z2...
            +(769/4346496000+(19182533/5339670336000000-8992463/1708694507520000000*z2)*z2...
            +(-280393/284782417920000-367176113/10252167045120000000*z2-7963/2531399270400000*z)*z)*z)*z)*z);
end
end
% plot de exacte en numerieke oplossingen ter illustratie van de
% (in)stabiliteit
%   @param max: bovengrens voor het interval waarin de oplossingen
%               vergeleken worden
%   @param inpute: de frequentie van de testvergelijking
function RKs3_instability_test(max,inpute)
global e
h = 0.125;
e = inpute/h;
% pas hiervoor theta aan in RKs3 (moet constant zijn)!
q = RKs3([0,max],h,1,1,@simple_1d_system,@(~,~) 8);
qe = exact_sol([0,max],h,'simple_1d_exact');
figure;
plot(0:h:max,q,0:h:max,qe);
legend('RKs3','Exacte oplossing');
end

% stelsel voor het eendimensionaal lineair autonoom probleem 
function y = simple_1d_system(q)
global e
y(1) = e*q(2);
end
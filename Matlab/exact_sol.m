% geeft de exacte oplossing van een probleem over een interval met een
% zekere stapgrootte
%   @param xspan: het interval waarin de oplossing bepaald dient te worden
%   @param h: de stapgrootte
%   @param f: de naam van het probleem
function y = exact_sol(xspan,h,f)
y = [feval(f,xspan(1)).'];
i = xspan(1);
while i < xspan(2)-h/2
    i = i+h;
    sol = feval(f,i).';
    y = [y,sol];
end
end

% exacte oplossing van lineair autonoom probleem met frequenties 1 en e
%   @param t: tijdsstap waarin de oplossing geëvalueerd wordt
function sol = simple_exact(t)
global e;
sol=[sin(t),sin(e*t),cos(t),e*cos(e*t)];
end

% exacte oplossing van het probleem van euler
%   @param t: tijdsstap waarin de oplossing geëvalueerd wordt
function sol = euler_exact(t)
[SN,CN,DN] = ellipj(t,0.51); 
sol=[sqrt(1.51)*SN,CN,DN];
end

% oplossing van het probleem van kepler, benaderd met machinenauwkeurigheid
%   @param t: tijdsstap waarin de oplossing geëvalueerd wordt
function sol = kepler_exact(t)
global e;
E=t;
del=1;
while abs(del)>10^(-10)
  del=-(t-E+e*sin(E))/(-1+e*cos(E));
  E=E+del;
end
c=cos(E);
s=sin(E);
noem=1/(1-e*c);
sq=sqrt(1-e^2);
sol=[c-e,sq*s,-s*noem,sq*noem*c];
end

% exacte oplossing van het geperturbeerd keplerprobleem
%   @param t: tijdsstap waarin de oplossing geëvalueerd wordt
function sol = pert_kepler_exact(t)
global e;
a=1+e;
c=cos(a*t);
s=sin(a*t);
sol=[c,s,-a*s,a*c];
end

% exacte oplossing van eendimensionaal lineair autonoom probleem
%   @param t: tijdsstap waarin de oplossing geëvalueerd wordt
function sol = simple_1d_exact(t)
global e;
sol=[exp(e*t)];
end
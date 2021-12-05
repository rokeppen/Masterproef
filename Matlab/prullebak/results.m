% accuracy options for solvers used in methods
global opt;
opt = optimoptions('fsolve','Display','Off','MaxIter',5000,'MaxFunEvals',5000000,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16);

% problems
simple_alphas(2,[0,10])
pert_kepler_alphas([0,10])
kepler_alphas([0,10])
euler_alphas([0,10])
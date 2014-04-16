function [y, v, grad, aSet, k, flag] = proximalGradient( x, fval, grad0, aSet0, mu, eta, stepTol, problem, verbose )

% Solves the U-step subproblem.
%
% For the bound constrained problem, this step is a projected
% gradient iteration:
%
%		y_{k+1} = Proj_Omega { y_k - alpha*grad_k },
%
% where Omega = { x : l \le x \le u }.
%

g = @problem.grad;
a = @problem.activeSet;

% Projected gradient step in the V-space.
v = fval;
grad = gUV(grad0,~aSet0);
a0 = (grad'*grad)/(grad'*problem.hprod(grad));
y = x;
k = 0;
maxIter = 15;
bestDecrease = -Inf;

exitGrad = false;
exitIter = false;

while ~exitGrad && ~ exitIter

vo = v;
k = k + 1;
    
[y, v, flag] = lineSearch(y, fval, grad0, -grad, mu, a0, stepTol, problem, verbose);

gr = g(y);
grad = gUV(gr,~aSet0);

decrease = vo-v;
bestDecrease = max(bestDecrease, decrease);
exitGrad = decrease <= eta*bestDecrease;
exitIter = k >= maxIter;

end

grad = gr;
aSet = a(y);

if exitGrad, flag = 'Sufficient decrease.'; end
if exitIter, flag = 'Max iterations.'; end

end % proximalGradient
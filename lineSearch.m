function [py, phik, flag, problem] = lineSearch( x, fval, g, d, mu, a0, stepTol, problem, verbose)

% Carries out a projected line search on x.

f = @problem.func;
p = @problem.prox;

alpha = a0;
phi0 = fval;
phik = phi0;

exitDesc = false;
exitStep = false;

pprint('\n', '', verbose);

while ~exitDesc && ~exitStep
    
    pprint('alpha = %f \n', alpha, verbose);
    
    y = x + alpha*d;
    py = p(y, alpha);
    phik = f(py);
    
    % Should be a descent step!
    if norm(py - y) < 1e-12 && g'*(py-x) >= 1e-12 
        fprintf('\n');
        error('Not taking a descent step.')
    end

    dxg = g'*(py-x); % directional derivative
    exitDesc = phik <= phi0 + mu*dxg; 
    exitStep = alpha < stepTol;
    
    alpha = 0.5*alpha;
    
    pprint('phik = %1.6f \n', phik, verbose);
    
end % while

if exitDesc, flag = 'Sufficient decrease.'; end
if exitStep, flag = 'Step size too small.'; py = x; phik = phi0; end

pprint('Reason for linesearch exit: %s \n', flag, verbose);

end % lineSearch
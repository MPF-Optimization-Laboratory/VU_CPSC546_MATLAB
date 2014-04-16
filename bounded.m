classdef bounded < Problem
    
    methods
        function prob = bounded(H, b, bl, bu)
            prob = prob@Problem(H, b);
            if isempty(bl)
                bl = zeros*ones(prob.n,1); % contrain to non-negative orthant
            end
            if isempty(bu)
                bu = Inf*ones(prob.n,1); % uncontrained above
            end
            if bl > bu
                fprintf('\n');
                error('Invalid lower and upper bounds.');
            end
            prob.bl = bl;
            prob.bu = bu;
            
        end
        
        %----------------------------------------------------------------
        function [f, prob] = func(prob, x)
          % Evaluates
          %                 q(x) = 0.5*x'Hx + b'x
          % at given (x).
          % Updates function eval counter and returns f = q(x).
          %--------------------------------------------------------------
            i = x - prob.bl < eps & prob.bu - x < eps;
            if sum(i) > 0
                indicator = inf;
            else
                indicator = 0;
            end
            Hx         = prob.hprod(x);
            f          = 0.5*(x'*Hx) + prob.b'*x + indicator;
        end % func
        
        %----------------------------------------------------------------
        function [g, prob] = grad(prob, x)
          % Evaluates 
          %                 q'(x) = Hx + b
          % at given (x).
          % Updates gradient eval counter and returns g = q'(x).
          %--------------------------------------------------------------
            Hx         = prob.hprod(x);
            g          = Hx + prob.b;
        end % grad
        
        %----------------------------------------------------------------
        function [x, prob] = prox(prob, x, ~)
          % Computes the projection mapping for problem.
          % Returns the mapped iteration.
          %--------------------------------------------------------------
            x = min(prob.bu, max(prob.bl,x));
        end % prox        
        
        %----------------------------------------------------------------
        
        function [aSet, prob] = activeSet(prob, x)
          % Calculates the active set of x.
          % Returns logical: union of lowAct (lower active) and
          %                           upAct  (upper active).
          %--------------------------------------------------------------
          lowAct = x == prob.bl;
          upAct  = x == prob.bu;
          aSet   = lowAct | upAct;
        end % activeSet
    end
end
        
        
classdef Problem
  % Problem: A class that defines a problem for vu.m
  %
  % Defines Problem class, includes a constructor function for problem.
  
    properties
        
      H              % Hessian
      b              % linear term
      bl             % lower bound
      bu             % upper bound
      n              % dimension
      explicitH      % NOTE: For now, assume H is given explicitly.
      
    end % properties.
    
    methods
      
        %----------------------------------------------------------------
        function prob = Problem(H, b)
          % Problem: Constructor for the problem class.
          % Creates a problem object for use in vu.m.
          %
          % Probelm form:
          %
          %     minimize   0.5*x'Hx + b'x
          %         s.t.   bl \le x \le bu         
          %--------------------------------------------------------------          
          prob.explicitH = true; % NOTE: For now, assume H is given explicitly.
          
          if prob.explicitH
              [N,m] = size(H);
              if N ~= m
                  fprintf('\n');
                  error('Given Hessian is not square.');
              end
              if ~all(all(tril(H)==triu(H).'))
                  fprintf('\n');
                  fprintf('Given Hessian is not symmetric.');
              end

              if isempty(b)
                  b = zeros(N,1);
              end
          end
              
          prob.H = opMatrix(H);
          prob.b = b;
          prob.n = N;

        end % problem.m
        
        %----------------------------------------------------------------
        function [z, prob] = hprod(prob, x)
          % Evaluates Hessian products with given (x).
          % Returns z = Hx.
          %--------------------------------------------------------------   
            
          if prob.explicitH
            z = prob.H*x;
          else
            z = prob.H(x,1);
          end

        end % hprod
 
        %----------------------------------------------------------------
        function [Z, ZHZ, prob] = reducedHess(prob, aSet)
          % Computes the reduced Hessian for solving reduced
          % subproblem in U-space:
          %--------------------------------------------------------------  
            
            Z   = Zop(aSet);
            ZHZ = ZHZop(prob.H, Z, aSet);
          
        end % reducedHess
  
    end % methods
    
end % Problem


function [ x, f, eflag, aSet, problem ] = vu( problem, xint, options )

% Implementation of VU-decomposition method for quadratic bound
% constrainted problems on the non-negative orthant:
%
%   [ xopt, fopt, eflag, aSet, problem ] = vu( problem, xint, options );
%
% Solves problems of the form
%
%               minimize   q(x) 
%                   s.t.   x in Omega,
%
% where
%    q(x)       is a (convex) quadratic, q(x) = 1\2 x'Hx + b'x, and  
%   Omega       is the non-negative orthant Omega = { x : lb \le x \le ub }.
%
% Input arguments:
%
%   problem     information for the given optimization problem.
%
%                   properties: H, b, Omega (lb, ub for bound constrained 
%                               problem), n, explicitH.
%
%                   methods:    func, grad, prox, activeSet, 
%                               hprod, reducedHess
%
%   xint        an estimate of the solution.
%
%   options     is a structure of options. Any unset options are
%               set to their default value; set options = [] to use all
%               default values.
%
% Output arguments:
%
%   x       is the optimal solution found by the algorithm
%   f       is the optimal function value found       
%   eflag   is the reason for the method terminating
%   aSet    is the optimal active set
%   problem is the problem structure with the terminating values

% IMPLEMENTED BY
% Julie A. Nutini
% University of British Columbia
% jnutini@cs.ubc.ca
%%----------------------------------------------------------------------
%% Set default input arguments.
%%----------------------------------------------------------------------
   if nargin < 1
      error('At least one argument is required.');
   end
   nogivenx = nargin < 2 || isempty(xint);
   if nargin < 3, options = []; end
%-----------------------------------------------------------------------
% Grab input options and set defaults where needed.
%-----------------------------------------------------------------------
  if isempty(options)
     options = as_setparms;
  end
  
  fid       = options.fid;
  callback  = isa(options.callback, 'function_handle');
  maxIterPG = options.maxIterPG;
  maxIterCG = options.maxIterCG;
  mu        = options.mu;
  eta        = options.eta;
  tau       = options.tau;
  relTol    = options.relTol;
  stepTol   = options.stepTol;
%-----------------------------------------------------------------------
% Initialize various variables.
%-----------------------------------------------------------------------
  itn = 0; % Iteration counter
  n = problem.n;
  nIterPG = 0;
  nIterCG = 0;
  verbose = 0; % Detailed printing.
  
% Exit flags.
  eflag            = 0;
  EXIT_OPTIMAL     = 1;
  EXIT_TOO_MANY_PG = 2;
  EXIT_TOO_MANY_CG = 3;
  EXIT_REQUESTED   = 4;
  exit_msg = {
      'OPTIMAL!'
      'Too many PG iterations.'
      'Too many CG iterations.'
      'iternal error'
      'Exit requested'
             };
%-----------------------------------------------------------------------
% Initialize cold/hot-start.
%-----------------------------------------------------------------------
  if nogivenx % Coldstart.
     x = zeros(n,1);
  else
     x = xint;
  end
  x = problem.prox(x, 1);
  g = problem.grad(x); % Used in optimality calculation.
%-----------------------------------------------------------------------
% Iteration log.
%-----------------------------------------------------------------------
  logHead = sprintf('\n %3s   %5s \t %4s \t\t %11s \t\t %6s \t%4s     %3s  %3s',...
            'Itn','PG/CG', 'Flag', 'f(x)','guNorm','|aSet|','nPG','nCG');
  logBody0 = '\n %2i \t %2s \t %20s \t  %-+01.4e \t %1.4e\t %6i    %4i  %4i';
  logBody1 = '\n %2i \t %2s \t %20s \t  %-+03.4e \t %1.4e\t %6i    %4i  %4i';
  fprintf(fid,'\n');
  fprintf(fid,logHead);
  fprintf(fid, ...
  '\n --------------------------------------------------------------------------------------------');
%-----------------------------------------------------------------------
% Begin iterations.
%-----------------------------------------------------------------------
  aSet = problem.activeSet(x);
  guNorm = norm(g);
  f = problem.func(x);
  fprintf(fid, logBody0, itn, '', '', f, guNorm, sum(aSet), nIterPG,nIterCG);
  
  while true
  %---------------------------------------------------------------------
  % Test exit conditions.
  %---------------------------------------------------------------------
    if nIterPG >= maxIterPG && eflag == 0
       eflag = EXIT_TOO_MANY_PG;
    end
    if nIterCG >= maxIterCG && eflag == 0
       eflag = EXIT_TOO_MANY_CG;
    end
    if callback
        reqexit = options.callback();
        if reqexit
            eflag = EXIT_REQUESTED;
        end
    end
    if guNorm <= tau
        eflag = EXIT_OPTIMAL;
    end

    % Act on exits
    if eflag ~= 0, break, end
    
    itn = itn + 1;
  %---------------------------------------------------------------------
  % Phase 1: Solve min_v f(x + 0 \oplus v) via proximal gradient.
  %---------------------------------------------------------------------   
    if sum(aSet) ~= 0 % else, default to CG as the problem is smooth
        [x, f, g, aSet, k, flagPG] = proximalGradient(x, f, g, aSet, mu, eta, stepTol, problem, verbose);
        nIterPG = nIterPG + k;
        gu = gUV(g,aSet);
        guNorm = norm(gu);
        fprintf(fid, logBody1, itn, 'PG', flagPG, f, guNorm, sum(aSet), nIterPG,nIterCG);
    end
    
  % Check optimality, max # of iterations.
    if guNorm <= tau || nIterPG >= maxIterPG
       continue
    end
  %---------------------------------------------------------------------
  % Phase 2: Explore U space via CG and solve: (H_u)u = -g_u.
  %---------------------------------------------------------------------  
  [Z, ZHZ] = problem.reducedHess(aSet);
  Ztg = Z(g,2);
  
  if sum(aSet) == 0 && itn ~= 1 % smooth problem
    [dZPCG, flagPCG, ~, CGiter, r] = pcg(ZHZ,-Ztg, tau, n);
  else % solve in U-space
     [dZPCG, flagPCG, ~, CGiter, r] = pcg(ZHZ,-Ztg, relTol, 30);
  end
  d = Z(dZPCG, 1);
  
  % Funny pcg quirk -- returns CGiter = 1 when flagPCG = 1.
  if flagPCG == 1
      nIterCG = nIterCG + length(r)-1;
  else
      nIterCG = nIterCG + CGiter;
  end
  
  [x, f, ~, problem] = lineSearch( x, f, g, d, mu, 1, stepTol, problem, verbose);
  g = problem.grad(x);
  aSet = problem.activeSet(x);
  gu = gUV(g,aSet);
  guNorm = norm(gu);

  fprintf(fid, logBody1, '', 'CG', pcgflag(flagPCG), f, guNorm, sum(aSet), nIterPG, nIterCG);
  end % while (outer iteration)
  
  %---------------------------------------------------------------------
  % Exit.
  %---------------------------------------------------------------------
  fprintf(fid, ...
  '\n --------------------------------------------------------------------------------------------');
  fprintf(fid, logBody1, '', '', 'FINISHED      ', f, guNorm, sum(aSet), nIterPG, nIterCG);
    fprintf(fid, ...
  '\n --------------------------------------------------------------------------------------------\n');
  fprintf(fid, ' %s \n', exit_msg{eflag});
  fprintf(fid, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % function driver_vu.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


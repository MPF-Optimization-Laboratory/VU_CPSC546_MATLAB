function [pf, mv_count, cond_n, fval_qp, fval_vu, time_qp, time_vu] = testProblems( prob, type )

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

% Type 1
M = [2800, 3600, 4400, 5200, 6000];
N = [2000, 2400, 2800, 3200, 3600];

% Type 2
sparsity = [0.0020, 0.0040, 0.0060, 0.0080, 0.0100, 0.0120];

if type == 1
    m = M(prob);
    n = N(prob);
    A = rand(m,n);
elseif type == 2
    m = 12000;
    n = 6400;
    A = sprand(m,n,sparsity(prob));
elseif type == 3
    load lp_pilot87
    A = Problem2.A;
    b1 = Problem2.b;
    n = size(A,2);
    xint = ones(n,1);
end

if type ~= 3
    b1 = rand(m,1);
    xint = zeros(n,1);
end

H = A'*A;
b = -A'*b1;
c = b1'*b1;
bl = zeros(n,1);
bu = Inf*ones(n,1);

if type == 1
    cond_n = cond(A);
else
    cond_n = 0;
end

% quadprog
tic;
[~,fval_qp] = quadprog(H,b,[],[],[],[],bl,bu,[],'TolFun',1e-6);
fval_qp = fval_qp + 0.5*c;
time_qp = toc;
% fval_qp = 0;
% time_qp = 0;

% Solve using VU heuristic
tic;
problem = bounded(H, b, bl, bu);
[~, fval_vu, eflag, ~] = vu(problem, xint);
fval_vu = fval_vu + 0.5*c;
time_vu = toc;
mv_count = problem.H.counter.mode1;

pf = eflag == 1 | eflag == 5; % eflag = 1 is optimal.
                              % eflag = 5 is not optimal. 

end


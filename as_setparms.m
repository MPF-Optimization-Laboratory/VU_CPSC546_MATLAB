function defopts = as_setparms(opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin default options.  Use only lower-case fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generic options
defopts.relTol     = 1e-6; % CG optimality tolerance 
defopts.callback   =    [];
defopts.fid        =     1;

% Options for gpcg submethods: conjGrad, gradProj, projSerach.
% Defaults taken from \cite{MT91}.
defopts.mu      =   0.3;	% ls sufficient decrease, eq. (4.1).
defopts.eta     =   0.25;	% ls sufficient decrease, eq. (4.1).
defopts.tau     =   1e-6;	% optimality tolerance (1e-4)
defopts.stepTol =   1e-10;	% Minimimu step length.

defopts.maxIterCG = 1000;
defopts.maxIterPG = 500;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End default options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick return if no user opts
if nargin == 0 || isempty(opts) 
    return
end
 
% List of valid field names
vfields = fieldnames( defopts );
 
% Grab valid fields from user's opts
for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
        defopts.(field) = opts.(field);
    end
end
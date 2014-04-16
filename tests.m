ttl = 'Test Problems';
T(:,1:2) = {
    'Test 1-1',@()testProblems(1,1)
    'Test 1-2',@()testProblems(2,1)
    'Test 1-3',@()testProblems(3,1)
    'Test 1-4',@()testProblems(4,1)
    'Test 1-5',@()testProblems(5,1)
    'Test 2-1',@()testProblems(1,2)
    'Test 2-2',@()testProblems(2,2)
    'Test 2-3',@()testProblems(3,2)
    'Test 2-4',@()testProblems(4,2)
    'Test 2-5',@()testProblems(5,2)
    'Test 2-6',@()testProblems(6,2)
%     'Test 3-1',@()testProblems(1,3)
    };

R = cell(size(T,1),8);
pf = {'**FAIL**','pass'};
for i=1:size(T,1);
    [passed, mv, cond_n, qp, vu, t_qp, t_vu] = T{i,2}();
    R(i,:) = {T{i,1}, pf{passed+1}, mv, cond_n qp, vu, t_qp, t_vu};
end
table(sprintf('%s - %i/%i',ttl,sum(strcmp(R(:,2),'pass')),size(R,1)),R,       ...
    'header',{'Test Name','pass/fail', '#MV prods (VU)', 'condition #', 'fval_qp','fval_vu', ...
    'quadprog (s)', 'VU (s)'},'format',{'%s','%s','%i', '%i', '%1.5e','%1.5e','%3.3g','%3.3g'});


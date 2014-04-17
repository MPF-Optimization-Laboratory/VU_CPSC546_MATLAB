VU_CPSC546_MATLAB
=================

MATLAB code implementing VU algorithm for convex quadratic function
restricted to the non-negative orthant.

The included files are described below:

Problem.m -- class, creates Problem structure
bounded.m -- subclass, defines specific methods for subclass of
             bounded convex quadratic problems.

vu.m               -- Implementation of the VU method.
as_setparms.m      -- File containing parameter values used.
proximalGradient.m -- Implementation of the projected gradient method.
lineSearch.m       -- Implementation of a projected linesearch.
ZHZop.m, Zop.m     -- Used in the construction of the reduced subproblem
                      in the U-space.
gUV.m              -- Forms the gradient in the U-space or the V-space, i.e.,
                      sets the dimensions not in the U-space or not in the 
                      V-space to 0.
pcgflag.m          -- Used to indicate stopping reason of pcg.m (MATLAB function).
pprint.m           -- Used with verbose = 1 input to print out detailed log.
table.m            -- Used to print out results table.


lp_pilot87.m    -- Real-world data set for future work. If the user uncomments 
                   'Test 3-1',@()testProblems(1,3) in tests.m, this problem
                   will run, but vu.m does not currentl converge.


TO RUN TEST PROBLEMS:
testProblems.m  -- List of test problems and how to form Problem and 
                   bounded objects, then run vu.m.
tests.m         -- RUN THIS to generate paper results.

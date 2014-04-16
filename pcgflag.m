function flag = pcgflag(flagPCG)
% Gets string description of pcg stopping flag.
% NOTE: flag == 2 corresponds to pre-conditioner.

    if flagPCG == 0
        flag = 'Optimal.';
    elseif flagPCG == 1
        flag = 'Max iterations.';
    elseif flagPCG == 3
        flag = 'Stagnated.';
    elseif flagPCG == 4
        flag = 'Scalar too small.';
    end

end
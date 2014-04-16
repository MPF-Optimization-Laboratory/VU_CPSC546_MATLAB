function Z = Zop(aSet)
% Computes 'reduced identity' according to aSet.

Z = @(d, mode)(Zop_internal(d, mode, aSet));

end

function y = Zop_internal(d, mode, aSet)
    n = length(aSet);
    if mode == 1
        y = zeros(n, 1);
        y(~aSet) = d;
    else
        y = d(~aSet);
    end
end

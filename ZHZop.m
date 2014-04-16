function ZHZ = ZHZop(H, Z, aSet)
% Computes reduced Hessian according to aSet.

ZHZ = @(d)(ZHZop_internal(d, H, Z, aSet));

end

function y = ZHZop_internal(d, H, Z, aSet)

    if sum(aSet) == 0
        y   = H*d;
    else
        Zd2  = Z(d, 1);
        HZd2 = H*Zd2;
        y   = Z(HZd2, 2);
    end

end
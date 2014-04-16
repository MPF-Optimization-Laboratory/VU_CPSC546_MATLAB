function gu = gUV(g,i)
% Takes g and sets either the U-space or V-space entires to 0.
% i = aSet (U-space) or i = ~aSet (V-space).

  gu = g;
  gu(i) = 0;
  
end
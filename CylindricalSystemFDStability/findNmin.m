function [val,ind] = findNmin(b,N)
[~,ind] = min(b);
 mval = max(b);
if N==0
    [val,ind] = min(b);
else
    bnew = b;
    bnew(ind) = mval;
    [val,ind] = findNmin(bnew,N-1);
end

end


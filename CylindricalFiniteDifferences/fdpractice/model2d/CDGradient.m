function [ ddf ] = CDGradient( dfhan,X0, h )

N = max(size(X0));
ddf = zeros(N);
%dfhan = @(X) sin(X);
parfor ik = 1:N
    loc = zeros(size(X0)); loc(ik) = h/2;
    dfhanloc = dfhan;
    %ddf(ik,:) = (dfhanloc(X0 + loc) - dfhanloc(X0 - loc))'/h;
    ddf(ik,:) = (dfhanloc(X0 + loc) - dfhanloc(X0 - loc))/h;

    
    
end
%ddf = sparse((1/2)*(ddf + ddf'));


end


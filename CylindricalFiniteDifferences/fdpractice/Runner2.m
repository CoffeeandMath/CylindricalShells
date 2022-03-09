Nnodes = 250;
Next = Nnodes + 1;
L = 1;
lvals = linspace(0,L,Nnodes);

h = lvals(2) - lvals(1);
lextvals = [lvals,lvals(end)+h];
ui = zeros(Next,1);


DD = tridiag(Next,-2,1,1)/h^2; DD(1,:) = 0;DD(1,1) = 1;
DDconst = DD(1:(end-1),:);

DDconststr = [DDconst;zeros(1,size(DDconst,2))];
DDconststr(end,end) = 1/h;
DDconststr(end,end-2) = -1/h;

Ffun = @(S) S;



tol = 10^-15;
err = 2*tol;
sc = 1e0;
cnt = 0;
Id = eye(size(DDconststr)); Id(1,1) = 0; Id(end,end) = 0;
fvalsext = [Ffun(lvals)';0];


L = DDconststr - Id;
usol = L\fvalsext;

usoliter = usol + 0.01*rand(size(usol));
usoliter = zeros(size(usoliter));
while err>tol
    cnt = cnt +1
    
    du = L\(L*usoliter-fvalsext);
    
   
    usoliter = usoliter - sc*du;
    err = norm(du)/max(size(du))
end




plot(lvals,usol(1:(end-1)),'.')
hold all
fanal = @(S) -(exp(1 - S) - exp(1 + S) + S + exp(2)*S)/(1 + exp(2));
plot(lvals,fanal(lvals))
plot(lextvals,usoliter)
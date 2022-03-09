Nnodes = 5000;
L = 1;
lvals = linspace(0,L,Nnodes);
h = lvals(2) - lvals(1);
ui = zeros(Nnodes,1);
dui = zeros(Nnodes,1);
ddui = zeros(Nnodes,1);
dddui = zeros(Nnodes,1);
ddddui = zeros(Nnodes,1);


D = tridiag(Nnodes,0,1,-1)/(2*h); D(1,1) = -1/h; D(1,2) = 1/h;
DD = tridiag(Nnodes,-2,1,1)/h^2; DD(1,1) = -1/h; D(1,2) = 1/h;
DDD = -D'*DD;
DDDD = DD'*DD;

Ffun = @(S) 0.5*sin(5*S);

dui = D*ui;
ddui = DD*ui;
dddui = DDD*ui;
ddddui = DDDD*ui;
ddfduduppi = zeros(Nnodes,1);
ddfdupduppi = zeros(Nnodes,1);
ddfdduppi = zeros(Nnodes,1);
dfdupi = zeros(Nnodes,1);
dfdui = zeros(Nnodes,1);

Fi = Ffun(lvals');
fprintf('constructed\n');
for i = 1:Nnodes
    ddfduduppi(i) = ddfdudupp(ui(i),dui(i),ddui(i));
    ddfdupduppi(i) = ddfdupdupp(ui(i),dui(i),ddui(i));
    ddfdduppi(i) = ddfddupp(ui(i),dui(i),ddui(i));
    dfdupi(i) = dfdup(ui(i),dui(i),ddui(i));
    dfdui(i) = dfdu(ui(i),dui(i),ddui(i));
end

L = D*(ddfduduppi.*dui + ddfdupduppi.*ddui) + (D*(ddfdduppi)).*dddui + ddfdduppi.*ddddui - D*dfdupi + dfdui - Fi;


tol = 10^-8;
err = 2*tol;
sc = 1e-2;
cnt = 0;
while err>tol
    cnt = cnt +1
    for i = 1:Nnodes
        ddfduduppi(i) = ddfdudupp(ui(i),dui(i),ddui(i));
        ddfdupduppi(i) = ddfdupdupp(ui(i),dui(i),ddui(i));
        ddfdduppi(i) = ddfddupp(ui(i),dui(i),ddui(i));
        dfdupi(i) = dfdup(ui(i),dui(i),ddui(i));
        dfdui(i) = dfdu(ui(i),dui(i),ddui(i));
    end
    
    L = D*(ddfduduppi.*dui + ddfdupduppi.*ddui) + (D*(ddfdduppi)).*dddui + ddfdduppi.*ddddui - D*dfdupi + dfdui - Fi;
    
    ui = ui - sc*L;
    err = norm(L)
end

plot(lvals,ui,'.')
hold all
plot(lvals,Fi)

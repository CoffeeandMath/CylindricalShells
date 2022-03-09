clear all
L = 1;

Nnodes = 50;
lvals = linspace(0,L,Nnodes);
h = lvals(2) - lvals(1);
lextvals = [-2* h, h, lvals, L + h, L + 2*h];

Next = Nnodes + 4;

Ndof = 2*Next;

uvec = zeros(Ndof,1);


ui = cell(Next,1);
for i = 1:max(size(ui))
    ui{i} = zeros(2,1);
    
    if (lextvals(i) > 0 && lextvals(i) < L)
        ui{i}(1) = lextvals(i)^7 * (L - lextvals(i))^7;
    else
        ui{i}(1) = 0;
    end
    
    
end

upi = calcup(ui,h,false);

uppi = calcupp(ui,h,false);

f = @(u,up,upp) 0.5*(norm(u)^2 + norm(up)^2 + norm(upp)^2);
dfdu = @(u,up,upp) u;
dfdup = @(u,up,upp) up;
dfdupp = @(u,up,upp) upp;

 
g2n = -1:(Nnodes+2);
n2g = @(S) S+2;


fi = zeros(Next,1);
dfdui = cell(Next,1);
dfdupi = cell(Next,1);
dfduppi = cell(Next,1);

for i = 2:(Next-1)
    fi(i) = f(ui{i},upi{i},uppi{i});
    dfdui{i} = dfdu(ui{i},upi{i},uppi{i});
    dfdupi{i} = dfdup(ui{i},upi{i},uppi{i});
    dfduppi{i} = dfdupp(ui{i},upi{i},uppi{i});
end

DfupDSi = calcup(dfdupi,h,true);
DDfuppDDSi = calcupp(dfduppi,h,true);

EL = [];

for i = 3:(max(size(DDfuppDDSi)) - 2)
     EL = [EL; DDfuppDDSi{i} - DfupDSi{i}];
end

ELExtend = [ui{1}; ui{2}; EL; ui{end-1}; ui{end}]
figure()
plot(ELExtend,'.')
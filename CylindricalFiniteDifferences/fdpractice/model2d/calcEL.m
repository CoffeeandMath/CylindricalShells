function ELExtend = calcEL(u)


L = 1;

Nnodes = size(u,1)/2 - 4;
lvals = linspace(0,L,Nnodes);
h = lvals(2) - lvals(1);
lextvals = [-2* h, h, lvals, L + h, L + 2*h];

Next = Nnodes + 4;

Ndof = 2*Next;

uvec = zeros(Ndof,1);


ui = setu(u);


upi = calcup(ui,h,false);

uppi = calcupp(ui,h,false);


 
g2n = -1:(Nnodes+2);
n2g = @(S) S+2;


fi = zeros(Next,1);
dfdui = cell(Next,1);
dfdupi = cell(Next,1);
dfduppi = cell(Next,1);

for i = 2:(Next-1)
    fi(i) = f(lextvals(i),ui{i},upi{i},uppi{i});
    dfdui{i} = dfdu(lextvals(i),ui{i},upi{i},uppi{i});
    dfdupi{i} = dfdup(lextvals(i),ui{i},upi{i},uppi{i});
    dfduppi{i} = dfdupp(lextvals(i),ui{i},upi{i},uppi{i});
end

DfupDSi = calcup(dfdupi,h,true);
DDfuppDDSi = calcupp(dfduppi,h,true);
DfuppDSi = calcup(dfduppi,h,true);

EL = [];


for i = 3:(max(size(DDfuppDDSi)) - 2)
     EL = [EL; DDfuppDDSi{i} - DfupDSi{i} + dfdui{i}];
end
EL = EL -0. * ones(size(EL));

EL(1:2) = ui{3} - [1;0];
EL((end-1):end) = ui{end-2} - [1;1];
%EL((end-1):end) = ui{end-2};
%ELExtend = [ui{1}; dfduppi{3}; EL; dfduppi{end-2}; dfdupi{end-2} - DfuppDSi{end-2}];
ELExtend = [ui{1}; upi{3}; EL; upi{end-2}; ui{end}];
end


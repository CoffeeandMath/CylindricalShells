clear all
L = 1;

Nnodes = 230;
lvals = linspace(0,L,Nnodes);
h = lvals(2) - lvals(1);
lextvals = [-2* h, -h, lvals, L + h, L + 2*h];

Next = Nnodes + 4;

Ndof = 2*Next;

uvec = zeros(Ndof,1);


ui = cell(Next,1);
for i = 1:max(size(ui))
    ui{i} = zeros(2,1);
    ui{i}(1) = 1;
    ui{i}(2) = lextvals(i);
    
    
    
end
ui{1} = zeros(2,1);
ufull = getu(ui);
Eltemp = calcEL(ufull);

tol = 1e-12;
err = 2*tol;
%sc = .000000002;
sc = 1;
cnt = 0;
while (err > tol)
    cnt = cnt + 1
    Eltemp = calcEL(ufull);
    DElhan = @(u) calcEL(u);
    %DElhannorm = @(u) norm(DElhan(u))^2;
    DEltemp = CDGradient(DElhan,ufull,1e-12);
    dL = DEltemp'\Eltemp;
    %dL = Eltemp;
    %dL = fminunc(DElhannorm,ufull);
    ufull = ufull - sc*dL;
    

    err = norm(dL)/size(dL,1)

    
end

[r,z] = getrz(ufull);


subplot(1,3,1)
plot(lvals',r(3:(end-2)),'.')
axis equal
hold all
subplot(1,3,2)
plot(lvals',z(3:(end-2)),'.')
axis equal
hold all
subplot(1,3,3)
hold all
plot(r(3:(end-2)),z(3:(end-2)),'b')
plot(-r(3:(end-2)),z(3:(end-2)),'b')
axis equal
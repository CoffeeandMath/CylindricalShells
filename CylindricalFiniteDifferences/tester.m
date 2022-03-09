close all
N = 2000;
nconst = 1;
L = 1.0;
lvals = linspace(0,L,N);
dl = lvals(2)-lvals(1);
Dm = zeros(N,N);


fval = 1.0;
Aval = 1.0;

% for i = 1:(N-1)
%     Dm(i,i+1) = 1/(2 * dl);
%     Dm(i+1,i) = -1/(2 * dl);
% end
% 
% Dm(1,1) = -1/dl;
% Dm(1,2) = 1/dl;
% 
% Dm(N,N) = 1/dl;
% Dm(N,N-1) = -1/dl;
% 
for i = 2:(N)
    Dm(i,i) = 1/dl;
    Dm(i,i-1) = -1/dl;
end
%Dm(N,N-1) = -1/dl;
Dm(N,N) = -1/dl;


Dm2 = (Dm')*Dm;

Dmv2 = zeros(N,N);

for i = 1:N
    Dmv2(i,i) = -2/(dl^2);
    if i < N
        Dmv2(i,i+1) = 1/(dl^2);
        Dmv2(i+1,i) = 1/(dl^2);
    end
end
Dmv2(N,N-1) = 1/dl;
Dmv2(N,N) = -1/dl;


[A,B] = eig(Dm2((nconst+1):end,(nconst+1):end));
figure()
hold all
for i = 1:2:10
    plot(lvals((nconst+1):end),A(:,i))
    amp = max(abs(A(:,i)));
    sgn = sign(A(2,i) - A(1,i));
    plot(lvals((nconst+1):end),amp*sgn*sin(i*pi*lvals((nconst+1):end)/2))
end


figure()

DDconst = Dm2((nconst+1):end,(nconst+1):end);
DDconst2 = Dmv2((nconst+1):end,(nconst+1):end);

fconst = fval*ones(N - nconst,1);

qsolve = (Aval*DDconst)\(fconst);
qsolve2 = (Aval*DDconst2)\(-fconst);
qext = [zeros(nconst,1);qsolve];
plot(lvals,qext,'b.')
hold all
%plot(lvals,[zeros(nconst,1);qsolve2],'r.')

C0 = fval*L/Aval;
fanal = @(S) -(fval/Aval)*S.^2 / 2 + C0 * S;
plot(lvals,fanal(lvals),'k-')

 dq = Dm*qext;

ddq = Dm2*qext;

figure()
subplot(1,2,1)
plot(lvals,dq)

hold all
subplot(1,2,2)
plot(lvals,ddq)
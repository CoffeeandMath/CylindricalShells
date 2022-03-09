N = 20;
Next = 2*N-1;
D = tridiag(Next,0,1,-1);

proj = zeros(N,Next);

for i = 1:N
    proj(i,2*i-1) = 1;
end

aveproj = zeros(N,Next);
for i = 1:N
    aveproj
end


uextend = zeros(Next,N);
for i = 1:N
    uextend(2*i-1,i) = 1;
    if i < N
        uextend(2*i,i) = 0.5;
        uextend(2*i,i+1) = 0.5;
    end
end
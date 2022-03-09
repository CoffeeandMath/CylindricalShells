function M = tridiag(N,a,b,c)
M = zeros(N,N);
M(  1:1+N:N*N) = a;
M(N+1:1+N:N*N) = b;
M(  2:1+N:N*N-N) = c;
end


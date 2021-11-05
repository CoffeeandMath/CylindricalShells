Ntsteps = 120;
Nnsteps = 5;

Ntvals = 1:Ntsteps;
Nnvals = 1:Nnsteps;

Neigs = 5;
eigsvals = zeros(Neigs,Nnsteps,Ntsteps);
for nt = 1:Ntsteps
    
    for nn = 1:Nnsteps
        file = ['build/matrices/ktildeout' , num2str(nn) , '_' , num2str(nt),'.csv']
        ktilde = readmatrix(file); 
        ktildesym = 0.5*(ktilde + ktilde'); 
        evs = eig(ktildesym);
        eigsvals(1:Neigs,nn,nt) = evs(1:Neigs);
    end
end
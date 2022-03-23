Nsteps = 800;
iter = 0:(Nsteps-1);
Neigs = 10;
lowesteigs = zeros(Neigs,Nsteps);


ploton = false;
eigvalplotoffset = false;
skip = 10;
%figure();
if not(ploton)
    skip = 1;
end
Neigskip = 0;
NFplot = 4;
noffset = 24;
Ndof = 6;
alleigs = zeros(Nsteps,NFplot+1,Neigs);


parfor i = 1:skip:Nsteps
    i
    tempmatrix = zeros(NFplot+1,Neigs);
    for fm = 0:NFplot
        Mv = load(['build/stabmatrices/stabmatrix' , num2str(iter(i)) , '_', num2str(fm) , '.csv']);
        [A,B] = eig(Mv);
        bv = diag(B);
        
        tempmatrix(fm+1,:) = bv(1:Neigs);
        
        subplot(1,NFplot+1,fm+1);
        
        
        
        
    end
    
    alleigs(i,:,:) = tempmatrix;
    
    

    
end


%%


figure();


for fm = 0:(NFplot)
    
    subplot(1,NFplot+1,fm+1)
    
    for ne = 1:Neigs
        plot(alleigs(:,fm+1,ne),'.')
        hold all
    end
    ylim([-0.1 2])
    xlim([1 625])
    title(num2str(fm))
    xlabel('Iteration')
    ylabel('Eigenvalue')
    
end



Nsteps = 800;
iter = 0:(Nsteps-1);
Neigs = 10;
lowesteigs = zeros(Neigs,Nsteps);


ploton = true;
padding = false;
padding2 = true;
eigvalplotoffset = false;
skip = 10;
figure();
if not(ploton)
    skip = 1;
end
Neigskip = 0;
Neigplot = 5;

Ndof = 2;
for i = 1:skip:Nsteps
    i
    Mv = load(['build/stabmatrices/stabmatrix' , num2str(iter(i)), '_0' , '.csv']);
    dl = max(size(Mv))/Ndof;
    [A,B] = eigs(Mv,Neigs,'SM');
    evs = sort(real(diag(B)));
    
    lowesteigs(:,i) = evs(1:Neigs);
    
    if ploton
        
        cntr = 0;
        for j = 0:(Neigplot-1)
            
            isviable = false;
            
            while (not(isviable))
                bv = real(diag(B));
                
                [bmin,bminind] = findNmin(bv,cntr);
                bmin
                evec = A(:,bminind);
                
                dr = diff(evec(1:2:end))/dl;
                dz = diff(evec(2:2:end))/dl;
                lnorm = max(abs(dr)+abs(dz))
                
                if (lnorm > 100)
                    cntr = cntr + 1
                    isviable = false;
                else
                    isviable = true;
                    cntr = cntr + 1
                end
                isviable = true;
            end
            
            rv = sortrows(readmatrix(['build/solutions/r_values_' num2str(iter(i)) '.csv']));
            zv = sortrows(readmatrix(['build/solutions/z_values_' num2str(iter(i)) '.csv']));
            subplot(Neigplot,2,2*j+1)
            
            sc = 1;
            hold off
            plot(rv(:,2),zv(:,2),'k.')
            hold all
            rpert = evec(1:2:end);
            zpert = evec(2:2:end);
            plot(rv(:,2) + sc*rpert ,zv(:,2) + sc*zpert,'r.');
            plot(rv(:,2) + sc*rpert ,zv(:,2) + sc*zpert,'r--');
            plot(rv(:,2) - sc*rpert,zv(:,2) - sc*zpert,'r.');
            plot(rv(:,2) - sc*rpert,zv(:,2) - sc*zpert,'r--');
            %         xlim([0 1])
            %         ylim([0 1])
            axis equal;
            title(num2str(bmin))
            subplot(Neigplot,2,2*j+2)
            
            
            offset = 1;
            hold off
            plot( real(evec(offset:2:(end-offset+1))),'r--')
            
            hold all
            plot(real( evec(offset:2:(end-offset+1))),'r.')
            
            plot( real(-evec(offset:2:(end-offset+1))),'r--')
            plot( real(-evec(offset:2:(end-offset+1))),'r.')
            plot( real(evec((offset+1):2:(end-offset+1))),'k--')
            plot(real( evec((offset+1):2:(end-offset+1))),'k.')
            plot( real(-evec((offset+1):2:(end-offset+1))),'k--')
            plot( real(-evec((offset+1):2:(end-offset+1))),'k.')
            xline(3);
            xline((size(evec,1))/2-2)
            xlim([0 size(evec,1)/2])
            
        end
        drawnow
    end
end

figure()
hold all
for i = 1:Neigs
    plot(iter,lowesteigs(i,:))
end

%% Plotting the deformation

ival = 0;

Mv = load(['build/stabmatrices/stabmatrix' , num2str(ival) , '.csv']);

[A,B] = eig(Mv);
bv = real(diag(B));
[bmin,bminind] = min(bv);
evec = A(:,bminind);

rv = sortrows(readmatrix(['build/solutions/r_values_' num2str(ival) '.csv']));
zv = sortrows(readmatrix(['build/solutions/z_values_' num2str(ival) '.csv']));

sc = 1.;
figure()
plot(rv(:,2),zv(:,2),'k.')
hold all
plot(rv(:,2)+ sc* evec(1:2:end),zv(:,2) + sc*evec(2:2:end));
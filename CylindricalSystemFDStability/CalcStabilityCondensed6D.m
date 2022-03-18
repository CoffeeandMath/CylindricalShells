Nsteps = 800;
iter = 0:(Nsteps-1);
Neigs = 10;
lowesteigs = zeros(Neigs,Nsteps);


ploton = false;
eigvalplotoffset = false;
skip = 10;
figure();
if not(ploton)
    skip = 1;
end
Neigskip = 0;
Neigplot = 4;
noffset = 12;
Ndof = 6;
parfor i = 1:skip:Nsteps
    i
    Mv = load(['build/stabmatrices/stabmatrix' , num2str(iter(i)) , '.csv']);
    dl = 1/size(Mv,1);
%     Mv = Mv((noffset+1):(end-noffset),(noffset+1):(end-noffset));
%     
%     Bv = eye(size(Mv));close
%     Bv(1:Ndof,1:Ndof) = 0;
%     Bv((end-Ndof+1):end,(end-Ndof+1):end) = 0;
%     [A,B] = eigs(Mv,Bv,Neigs,'SM');
    [A,B] = eig(Mv);
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
            plot(rv(2:(end-1),2) + sc*rpert ,zv(2:(end-1),2) + sc*zpert,'r.');
            plot(rv(2:(end-1),2) + sc*rpert ,zv(2:(end-1),2) + sc*zpert,'r--');
            plot(rv(2:(end-1),2) - sc*rpert,zv(2:(end-1),2) - sc*zpert,'r.');
            plot(rv(2:(end-1),2) - sc*rpert,zv(2:(end-1),2) - sc*zpert,'r--');
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
    plot(iter(1:skip:end),lowesteigs(i,1:skip:end))
end

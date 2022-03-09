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
Neigplot = 8;

Ndof = 6;
for i = 1:skip:Nsteps
    i
    Mv = load(['stabmatrices/stabmatrix' , num2str(iter(i)) , '.csv']);

%     Mv(3:4,:) = 0; 
%     O3 = [-0.5*eye(2),eye(2),zeros(2),-eye(2),0.5*eye(2)]; Mv(3:4,1:10) = O3;
%     Mv((end-3):(end-2),:) = 0; Mv((end-3):(end-2),(end-9):(end)) = O3;  
    dl = 1/size(Mv,1);
    Bv = eye(size(Mv)); Bv(1:2*Ndof,1:2*Ndof)= 0; Bv((end-2*Ndof+1):end,(end-2*Ndof+1):end) = 0;

    
    
    
    %         Mv(5,:) = 0;
    %         Mv(5,5) = 1;
    %
    %         Mv(7,:) = 0;
    %         Mv(7,7) = 1;
    %         Mv(8,:) = 0;
    %         Mv(8,8) = 1;
    [A,B] = eig(Mv,Bv);
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
            
            rv = sortrows(readmatrix(['solutions/r_values_' num2str(iter(i)) '.csv']));
            zv = sortrows(readmatrix(['solutions/z_values_' num2str(iter(i)) '.csv']));
            subplot(Neigplot,2,2*j+1)
            
            sc = .1;
            hold off
            plot(rv(:,2),zv(:,2),'k.')
            hold all
            rpert = evec(5:2:(end-4));
            zpert = evec(6:2:(end-4));
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
for i = Neigskip:Neigs
    plot(iter,lowesteigs(i,:))
end

%% Plotting the deformation

ival = 0;

Mv = load(['stabmatrices/stabmatrix' , num2str(ival) , '.csv']);

[A,B] = eig(Mv);
bv = real(diag(B));
[bmin,bminind] = min(bv);
evec = A(:,bminind);

rv = sortrows(readmatrix(['solutions/r_values_' num2str(ival) '.csv']));
zv = sortrows(readmatrix(['solutions/z_values_' num2str(ival) '.csv']));

sc = 1.;
figure()
plot(rv(:,2),zv(:,2),'k.')
hold all
plot(rv(:,2)+ sc* evec(1:2:end),zv(:,2) + sc*evec(2:2:end));
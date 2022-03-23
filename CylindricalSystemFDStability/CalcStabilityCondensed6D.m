Nsteps = 145;
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
noffset = 24;
Ndof = 6;
for i = 1:skip:Nsteps
    i
    Mv = load(['build/stabmatrices/stabmatrix' , num2str(iter(i)) , '.csv']);
    dl = 1/size(Mv,1);
    Mv = Mv((noffset+1):(end-noffset),(noffset+1):(end-noffset));
    Mv(2,2) = 1;
    Mv(5,5) = 1;
    m1 = Mv(1:6,1:6); m1 = eye(size(m1));
    m2 = Mv(1:6,7:end);
    m3 = Mv(7:end,1:6);
    Msub1 = Mv(7:end,7:end) - m3*(m1\m2);
    m1 = Msub1((end-5):end,(end-5):end);m1 = eye(size(m1));
    m2 = Msub1((end-5):end,1:(end-6));
    m3 = Msub1(1:(end-6),(end-5):end);
    Msub1 = 0.5*(Msub1 + Msub1');

    Msub2 = Msub1(1:(end-6),1:(end-6)) - m3*(m1\m2);
    Msub2 = 0.5*(Msub2+Msub2');

        Bv = eye(size(Mv));
        Bv(1:Ndof,1:Ndof) = 0;
        Bv((end-Ndof+1):end,(end-Ndof+1):end) = 0;
        [A,B] = eigs(Mv,Bv,Neigs,'SM');
    %[A,B] = eig(Mv);
%%
    ikk = 1;
    figure()
    subplot(1,2,1)
    plot(A(1:6:end,ikk))
    hold on
    plot(A(2:6:end,ikk))
    plot(A(3:6:end,ikk))
    legend('v_r','v_\theta','v_z')
    subplot(1,2,2)
    plot(A(4:6:end,ikk))
    hold on
    plot(A(5:6:end,ikk))
    plot(A(6:6:end,ikk))
    legend('w_r','w_\theta','w_z')
%%


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

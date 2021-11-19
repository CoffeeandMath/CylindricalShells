%% Plotting Symmetric Stability Results
close all;
figure();
evall = readmatrix('CylindricalSystem/build/evals.csv');
normalized = 0;
nfact = 1;
hold all
for i = 1:size(evall,2)
    if normalized
        nfact = abs(evall(1,i));
    end
    plot(evall(1:end,i)/nfact);
end
ylim([-1,1]*10^-5)
%% Plotting Asymmetric Stability Results

figure();
nfmode = 5;



for i = 1:nfmode
    subplot(1,5,i)
    hold all
    ev = readmatrix(['Cylindrical3DStability/build/evals' num2str(i) '.csv']);
    ev(abs(ev)<10^-10) = NaN;
    jstart = 1;
    if (i==1)
        jstart = 1;
    end
    for j = jstart:size(ev,2)
        if normalized
            nfact = abs(ev(1,j));
        else
            nfact = 1;
        end
        if mean(abs(ev(:,j))/size(ev(:,j),1) > 10^-10)
            plot(ev(:,j)/nfact,'k.')
            %ylim([-1,2])
            %ylim([-1,1]*10^-3)
        else
            fprintf('skipped \n')
        end
       
        
    end
    title(i)
end
subplot(1,5,1)
%ylim([-1,1]*10^-2)


%% Plotting Symmetric Stability Results
close all;
figure();
evall = readmatrix('CylindricalSystem/build/evals.csv');

hold all
for i = 1:size(evall,2)
    plot(evall(:,i)/evall(1,i));
end

%% Plotting Asymmetric Stability Results

figure();
nevals = 5;
normalized = 1;
nfact = 1;
for i = 1:nevals
    subplot(1,5,i)
    hold all
    ev = readmatrix(['Cylindrical3DStability/build/evals' num2str(i) '.csv']);
    jstart = 1;
    if (i==1)
        jstart = 1;
    end
    for j = 3:size(ev,2)
        if normalized
            nfact = abs(ev(1,j));
        else 
            nfact = 1;
        end
        plot(ev(:,j)/nfact)
        ylim([0,2])
    end
    title(i)
end
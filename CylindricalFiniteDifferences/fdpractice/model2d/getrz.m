function [r,z] = getrz(ufull)

    ui = setu(ufull);
    r = zeros(size(ui));
    z = zeros(size(ui));
    for i = 1:size(ui,1)
        r(i) = ui{i}(1);
        z(i) = ui{i}(2);
    end

end


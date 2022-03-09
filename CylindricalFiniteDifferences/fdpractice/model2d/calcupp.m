function uppi = calcupp(ui,h,inner)
    uppi = cell(size(ui));
    if inner
        range = 3:(max(size(ui))-2);
    else
        range = 2:(max(size(ui))-1);
    end
    for i = range
        uppi{i} = (ui{i+1} - 2 * ui{i} + ui{i-1})/(h^2);
    end

end

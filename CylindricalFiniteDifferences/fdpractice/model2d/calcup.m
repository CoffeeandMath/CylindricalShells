function upi = calcup(ui,h,inner)
    upi = cell(size(ui));
    if inner
        range = 3:(max(size(ui))-2);
    else
        range = 2:(max(size(ui))-1);
    end
    for i = range
        upi{i} = (ui{i+1} - ui{i-1})/(2*h);
    end

end

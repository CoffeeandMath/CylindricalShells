function ui = setu(ufull)
Ncell = size(ufull,1)/2;
ui = cell(Ncell,1);
cnt = 1;
for i = 1:Ncell
    ui{i} = zeros(2,1);
    ui{i}(1) = ufull(cnt);
    cnt = cnt + 1;
    ui{i}(2) = ufull(cnt);
    cnt = cnt + 1;
end


end


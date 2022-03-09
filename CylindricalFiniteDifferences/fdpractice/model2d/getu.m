function [ufull] = getu(ui)
ufull = zeros(2*size(ui,1),1);

cnt = 1;
for i = 1:size(ui,1)
    ufull(cnt) = ui{i}(1);
    cnt = cnt + 1;
    ufull(cnt) = ui{i}(2);
    cnt = cnt + 1;
end

end


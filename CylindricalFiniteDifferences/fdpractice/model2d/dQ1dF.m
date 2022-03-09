function val = dQ1dF(F,E,nu)

%val = 0.5*E*(tr(F)^2 - 2*(1 - nu)*det(F));
val = zeros(size(F));
val(1,1) = E * (F(1,1) + F(2,2) * nu);
val(1,2) = E * F(2,1) * (1 - nu);
val(2,1) = E * F(1,2) * (1 - nu);
val(2,2) = E * (F(2,2) + F(1,1) * nu);
end


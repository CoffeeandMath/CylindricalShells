function val = Q1(F,E,nu)

val = 0.5*E*(trace(F)^2 - 2*(1 - nu)*det(F));

end


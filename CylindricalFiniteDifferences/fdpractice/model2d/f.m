function val = f(S,u,up,upp)

% h = 0.1;
% Ev = 1e5;
% nu = 0.4;
% hsc = h^2/(12*(1 - nu^2));
% 
% 
% R = 1;
% 
% Ac = [1 , 0; 0, R^2];
% AC = [1, 0; 0, 1/R^2];
% 
% P = zeros(2,2); 
% P(1,1) = sqrt(Ac(1,1));
% P(1,2) = Ac(1,2)/sqrt(Ac(1,1));
% P(2,2) = sqrt(Ac(2,2) - Ac(1,2)^2/Ac(1,1));
% 
% ac = zeros(2,2);
% 
% a1 = [up(1);0;up(2)];
% a2 = [0;u(1);0];
% ac(1,1) = a1'*a1;
% ac(1,2) = a1'*a2;
% ac(2,1) = ac(1,2);
% ac(2,2) = a2'*a2;
% 
% aC = inv(ac);
% 
% E = 0.5*P*(AC*ac*AC - AC)*(P');
% 
% 
% 
% fv = cross(a1,a2);
% n = fv/norm(fv);
% 
% bc = zeros(2,2);
% ddxddS = [upp(1);0;upp(2)];
% ddxdSdtheta = [0;up(1);0];
% ddxddtheta = [-u(1);0;0];
% 
% bc(1,1) = n'*ddxddS;
% bc(1,2) = n'*ddxdSdtheta;
% bc(2,1) = bc(1,2);
% bc(2,2) = n'*ddxddtheta;
% 
% 
% B = P*aC*bc*aC*(P');
% 
% 
% 
% 
% val = Q1(E,Ev,nu) + hsc*Q1(B - [0,0;0,1],Ev,nu);
val = 0.5*(20 * norm(u)^2 + norm(up)^2 + norm(upp)^2);
end


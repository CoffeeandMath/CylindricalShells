function val = dfdupp(S,u,up,upp)


% 
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
% val = zeros(size(u));
% 
% DQ1 = dQ1dF(E,Ev,nu);
% DQ2 = dQ1dF(B - [0,0;0,1],Ev,nu);
% 
% for i = 1:size(u,1)
%     diupp = zeros(size(u)); diupp(i) = 1;
%     
%     dia1 = zeros(3,1);
%     dia2 = zeros(3,1);
%     
%     diac = zeros(2,2);
%     diac(1,1) = 2*dia1'*a1;
%     diac(1,2) = dia1'*a2 + a1'*dia2;
%     diac(2,1) = diac(1,2);
%     diac(2,2) = 2*dia2'*a2;
%     
%     diaC = - aC*diac*aC;
%     
%     diE = 0.5*P*(AC*diac*AC)*(P');
% 
%     
%     DiddxddS = [diupp(1);0;diupp(2)];
%     DiddxdSdtheta = zeros(3,1);
%     Diddxddtheta = zeros(3,1);
%     
%     difv = cross(dia1,a2) + cross(a1,dia2);
%     
%     din = (eye(3) - n*n')*difv / norm(fv);
%     
%     dib1 = zeros(2,2);
%     dib1(1,1) = din'*ddxddS;
%     dib1(1,2) = din'*ddxdSdtheta;
%     dib1(2,1) = dib1(1,2);
%     dib1(2,2) = din'*ddxddtheta;
%     
%     dib2 = zeros(2,2);
%     dib2(1,1) = n'*DiddxddS;
%     dib2(1,2) = n'*DiddxdSdtheta;
%     dib2(2,1) = dib2(1,2);
%     dib2(2,2) = n'*Diddxddtheta;
%     
%     
%     
%         
%     dibc = dib1 + dib2;
%     
%     diB = P*(aC*dibc*aC + 2*symmetric(aC*bc*diaC))*(P');
%     
%     val(i) = trace(DQ1*diE) + hsc*trace(DQ2*diB);
%     
% end

val = upp;

end


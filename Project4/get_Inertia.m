function Iner = get_Inertia(th2)

global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tq1 tq2;


M11 = Iz1 + Iz2 - L1^2*m1 + L1^2*m2 - L2^2*m2 + 2*L1*m1*r1 + 2*L2*m2*r2 + 2*L1*m2*r2*cos(th2);
 
 
M12 = - m2*L2^2 + 2*m2*r2*L2 + Iz2 + L1*m2*r2*cos(th2);
 
 
M21 = - m2*L2^2 + 2*m2*r2*L2 + Iz2 + L1*m2*r2*cos(th2);
 
 
M22 = - m2*L2^2 + 2*m2*r2*L2 + Iz2;



Iner = [M11, M12; M21, M22];

end


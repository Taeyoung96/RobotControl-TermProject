function Iner = get_Inertia(th2, th3)

global Iz1 Iz2 Iz3 L1 L2 L3 m1 m2 m3 r1 r2 r3;


M11 = Iz1 + Iz2 + Iz3 - L1^2*m1 + L1^2*m2 + L1^2*m3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L1*m1*r1 + 2*L2*m2*r2 + 2*L3*m3*r3 + 2*L1*m3*r3*cos(th2 + th3) + 2*L1*L2*m3*cos(th2) + 2*L1*m2*r2*cos(th2) + 2*L2*m3*r3*cos(th3);
M12 = Iz2 + Iz3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L2*m2*r2 + 2*L3*m3*r3 + L1*m3*r3*cos(th2 + th3) + L1*L2*m3*cos(th2) + L1*m2*r2*cos(th2) + 2*L2*m3*r3*cos(th3);
M13 = - m3*L3^2 + 2*m3*r3*L3 + Iz3 + L1*m3*r3*cos(th2 + th3) + L2*m3*r3*cos(th3);

M21 = Iz2 + Iz3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L2*m2*r2 + 2*L3*m3*r3 + L1*m3*r3*cos(th2 + th3) + L1*L2*m3*cos(th2) + L1*m2*r2*cos(th2) + 2*L2*m3*r3*cos(th3);
M22 =  Iz2 + Iz3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L2*m2*r2 + 2*L3*m3*r3 + 2*L2*m3*r3*cos(th3);
M23 = - m3*L3^2 + 2*m3*r3*L3 + Iz3 + L2*m3*r3*cos(th3);

M31 = - m3*L3^2 + 2*m3*r3*L3 + Iz3 + L1*m3*r3*cos(th2 + th3) + L2*m3*r3*cos(th3);
M32 =   - m3*L3^2 + 2*m3*r3*L3 + Iz3 + L2*m3*r3*cos(th3);
M33 =  - m3*L3^2 + 2*m3*r3*L3 + Iz3;


Iner = [M11, M12, M13; M21, M22, M23; M31, M32, M33;];

end
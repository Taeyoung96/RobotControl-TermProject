function get_G = get_Gravity(th1,th2)

global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tq1 tq2;


get_G = [g*(m1*r1*cos(th1) + m2*r2*cos(th1 + th2) + L1*m2*cos(th1));
                                    g*m2*r2*cos(th1 + th2)];

end


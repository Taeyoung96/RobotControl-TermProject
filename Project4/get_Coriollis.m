function Get_C = get_Coriollis(th2,dth1,dth2)


 global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tq1 tq2;


Get_C = [-L1*dth2*m2*r2*sin(th2)*(2*dth1 + dth2);
               L1*dth1^2*m2*r2*sin(th2)];

end


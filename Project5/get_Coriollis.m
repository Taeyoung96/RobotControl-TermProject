function Get_C = get_Coriollis(th2, th3, dth1,dth2,dth3)

 global L1 L2 L3 m1 m2 m3 r1 r2 r3 tq1 tq2;
                                                                                                                                                                                                                                                                 
C11 = - L1*dth2^2*m3*r3*sin(th2 + th3) - L1*dth3^2*m3*r3*sin(th2 + th3) - L1*L2*dth2^2*m3*sin(th2) - L1*dth2^2*m2*r2*sin(th2) - L2*dth3^2*m3*r3*sin(th3) - 2*L1*dth1*dth2*m3*r3*sin(th2 + th3) - 2*L1*dth1*dth3*m3*r3*sin(th2 + th3) - 2*L1*dth2*dth3*m3*r3*sin(th2 + th3) - 2*L1*L2*dth1*dth2*m3*sin(th2) - 2*L1*dth1*dth2*m2*r2*sin(th2) - 2*L2*dth1*dth3*m3*r3*sin(th3) - 2*L2*dth2*dth3*m3*r3*sin(th3);
C21 = L1*dth1^2*m3*r3*sin(th2 + th3) + L1*L2*dth1^2*m3*sin(th2) + L1*dth1^2*m2*r2*sin(th2) - L2*dth3^2*m3*r3*sin(th3) - 2*L2*dth1*dth3*m3*r3*sin(th3) - 2*L2*dth2*dth3*m3*r3*sin(th3);
C31 =  m3*r3*(L1*dth1^2*sin(th2 + th3) + L2*dth1^2*sin(th3) + L2*dth2^2*sin(th3) + 2*L2*dth1*dth2*sin(th3));

Get_C = [C11; C21; C31];

end

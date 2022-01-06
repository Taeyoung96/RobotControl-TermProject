function dydt = three_links(t,y)

global Iz1 Iz2 Iz3 L1 L2 L3 dth1 dth2 dth3 g m1 m2 m3 r1 r2 r3 tq1 tq2 tq3 th1 th2 th3

th1 = y(1);
dth1 = y(2);
th2 = y(3);
dth2 = y(4);
th3 = y(5);
dth3 = y(6);


dydt = [dth1;((-tq2+L2.*g.*m3.*cos(th1+th2)+g.*m2.*r2.*cos(th1+th2)+g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth1.^2.*m3.*r3.*sin(th2+th3)+L1.*L2.*dth1.^2.*m3.*sin(th2)+L1.*dth1.^2.*m2.*r2.*sin(th2)-L2.*dth3.^2.*m3.*r3.*sin(th3)-L2.*dth1.*dth3.*m3.*r3.*sin(th3).*2.0-L2.*dth2.*dth3.*m3.*r3.*sin(th3).*2.0).*(Iz2.*Iz3-L2.^2.*L3.^2.*m3.^2-Iz3.*L2.^2.*m2-Iz2.*L3.^2.*m3+Iz3.*L2.^2.*m3+L2.^2.*L3.^2.*m2.*m3+L2.^2.*L3.*m3.^2.*r3.*2.0+Iz3.*L2.*m2.*r2.*2.0+Iz2.*L3.*m3.*r3.*2.0-L2.^2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.*L3.^2.*m3.^2.*cos(th2)+Iz3.*L1.*L2.*m3.*cos(th2)-L2.*L3.^2.*m2.*m3.*r2.*2.0-L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz3.*L1.*m2.*r2.*cos(th2)-L1.*L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3)+L1.*L2.*L3.*m3.^2.*r3.*cos(th2).*2.0-L1.*L3.^2.*m2.*m3.*r2.*cos(th2)+L2.*L3.*m2.*m3.*r2.*r3.*4.0+L1.*L3.*m2.*m3.*r2.*r3.*cos(th2).*2.0))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)+((Iz2.*Iz3-L2.^2.*L3.^2.*m3.^2-Iz3.*L2.^2.*m2-Iz2.*L3.^2.*m3+Iz3.*L2.^2.*m3+L2.^2.*L3.^2.*m2.*m3+L2.^2.*L3.*m3.^2.*r3.*2.0+Iz3.*L2.*m2.*r2.*2.0+Iz2.*L3.*m3.*r3.*2.0-L2.^2.*m3.^2.*r3.^2.*cos(th3).^2-L2.*L3.^2.*m2.*m3.*r2.*2.0-L2.^2.*L3.*m2.*m3.*r3.*2.0+L2.*L3.*m2.*m3.*r2.*r3.*4.0).*(tq1-L2.*g.*m3.*cos(th1+th2)-g.*m2.*r2.*cos(th1+th2)-L1.*g.*m2.*cos(th1)-L1.*g.*m3.*cos(th1)-g.*m1.*r1.*cos(th1)-g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth2.^2.*m3.*r3.*sin(th2+th3)+L1.*dth3.^2.*m3.*r3.*sin(th2+th3)+L1.*L2.*dth2.^2.*m3.*sin(th2)+L1.*dth2.^2.*m2.*r2.*sin(th2)+L2.*dth3.^2.*m3.*r3.*sin(th3)+L1.*dth1.*dth2.*m3.*r3.*sin(th2+th3).*2.0+L1.*dth1.*dth3.*m3.*r3.*sin(th2+th3).*2.0+L1.*dth2.*dth3.*m3.*r3.*sin(th2+th3).*2.0+L1.*L2.*dth1.*dth2.*m3.*sin(th2).*2.0+L1.*dth1.*dth2.*m2.*r2.*sin(th2).*2.0+L2.*dth1.*dth3.*m3.*r3.*sin(th3).*2.0+L2.*dth2.*dth3.*m3.*r3.*sin(th3).*2.0))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)-(L1.*(-tq3+g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth1.^2.*m3.*r3.*sin(th2+th3)+L2.*dth1.^2.*m3.*r3.*sin(th3)+L2.*dth2.^2.*m3.*r3.*sin(th3)+L2.*dth1.*dth2.*m3.*r3.*sin(th3).*2.0).*(-Iz2.*m3.*r3.*cos(th2+th3)+Iz3.*L2.*m3.*cos(th2)+Iz3.*m2.*r2.*cos(th2)-L2.^2.*m3.^2.*r3.*cos(th2+th3)-L2.*L3.^2.*m3.^2.*cos(th2)+L2.^2.*m2.*m3.*r3.*cos(th2+th3)+L2.*L3.*m3.^2.*r3.*cos(th2).*2.0-L3.^2.*m2.*m3.*r2.*cos(th2)-L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3)+L2.^2.*m3.^2.*r3.*cos(th2).*cos(th3)-L2.*m2.*m3.*r2.*r3.*cos(th2+th3).*2.0+L3.*m2.*m3.*r2.*r3.*cos(th2).*2.0+L2.*m2.*m3.*r2.*r3.*cos(th2).*cos(th3)))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0);dth2;((-tq3+g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth1.^2.*m3.*r3.*sin(th2+th3)+L2.*dth1.^2.*m3.*r3.*sin(th3)+L2.*dth2.^2.*m3.*r3.*sin(th3)+L2.*dth1.*dth2.*m3.*r3.*sin(th3).*2.0).*(Iz1.*Iz3-L1.^2.*L3.^2.*m3.^2-Iz3.*L1.^2.*m1+Iz3.*L1.^2.*m2-Iz1.*L3.^2.*m3+Iz3.*L1.^2.*m3+L1.^2.*L3.^2.*m1.*m3-L1.^2.*L3.^2.*m2.*m3+L1.^2.*L3.*m3.^2.*r3.*2.0-L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2+Iz3.*L1.*m1.*r1.*2.0+Iz1.*L3.*m3.*r3.*2.0-L1.*L2.^2.*m3.^2.*r3.*cos(th2+th3)-L1.*L2.*L3.^2.*m3.^2.*cos(th2)+L1.^2.*L2.*m3.^2.*r3.*cos(th3)-Iz2.*L1.*m3.*r3.*cos(th2+th3)+Iz3.*L1.*L2.*m3.*cos(th2)-L1.*L3.^2.*m1.*m3.*r1.*2.0-L1.^2.*L3.*m1.*m3.*r3.*2.0+L1.^2.*L3.*m2.*m3.*r3.*2.0+Iz3.*L1.*m2.*r2.*cos(th2)+Iz1.*L2.*m3.*r3.*cos(th3)-L1.^2.*L2.*m3.^2.*r3.*cos(th2+th3).*cos(th2)-L1.*L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3)+L1.*L2.^2.*m3.^2.*r3.*cos(th2).*cos(th3)+L1.*L2.^2.*m2.*m3.*r3.*cos(th2+th3)+L1.*L2.*L3.*m3.^2.*r3.*cos(th2).*2.0-L1.*L3.^2.*m2.*m3.*r2.*cos(th2)-L1.^2.*L2.*m1.*m3.*r3.*cos(th3)+L1.^2.*L2.*m2.*m3.*r3.*cos(th3)+L1.*L3.*m1.*m3.*r1.*r3.*4.0-L1.^2.*m2.*m3.*r2.*r3.*cos(th2+th3).*cos(th2)-L1.*L2.*m2.*m3.*r2.*r3.*cos(th2+th3).*2.0+L1.*L2.*m1.*m3.*r1.*r3.*cos(th3).*2.0+L1.*L3.*m2.*m3.*r2.*r3.*cos(th2).*2.0+L1.*L2.*m2.*m3.*r2.*r3.*cos(th2).*cos(th3)))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)-((-tq2+L2.*g.*m3.*cos(th1+th2)+g.*m2.*r2.*cos(th1+th2)+g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth1.^2.*m3.*r3.*sin(th2+th3)+L1.*L2.*dth1.^2.*m3.*sin(th2)+L1.*dth1.^2.*m2.*r2.*sin(th2)-L2.*dth3.^2.*m3.*r3.*sin(th3)-L2.*dth1.*dth3.*m3.*r3.*sin(th3).*2.0-L2.*dth2.*dth3.*m3.*r3.*sin(th3).*2.0).*(Iz1.*Iz3+Iz2.*Iz3-L1.^2.*L3.^2.*m3.^2-L2.^2.*L3.^2.*m3.^2-Iz3.*L1.^2.*m1+Iz3.*L1.^2.*m2-Iz1.*L3.^2.*m3+Iz3.*L1.^2.*m3-Iz3.*L2.^2.*m2-Iz2.*L3.^2.*m3+Iz3.*L2.^2.*m3+L1.^2.*L3.^2.*m1.*m3-L1.^2.*L3.^2.*m2.*m3+L2.^2.*L3.^2.*m2.*m3+L1.^2.*L3.*m3.^2.*r3.*2.0+L2.^2.*L3.*m3.^2.*r3.*2.0-L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2+Iz3.*L1.*m1.*r1.*2.0+Iz3.*L2.*m2.*r2.*2.0+Iz1.*L3.*m3.*r3.*2.0+Iz2.*L3.*m3.*r3.*2.0-L2.^2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.*L3.^2.*m3.^2.*cos(th2).*2.0+Iz3.*L1.*L2.*m3.*cos(th2).*2.0-L1.*L3.^2.*m1.*m3.*r1.*2.0-L1.^2.*L3.*m1.*m3.*r3.*2.0-L2.*L3.^2.*m2.*m3.*r2.*2.0+L1.^2.*L3.*m2.*m3.*r3.*2.0-L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz3.*L1.*m2.*r2.*cos(th2).*2.0-L1.*L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3).*2.0+L1.*L2.*L3.*m3.^2.*r3.*cos(th2).*4.0-L1.*L3.^2.*m2.*m3.*r2.*cos(th2).*2.0+L1.*L3.*m1.*m3.*r1.*r3.*4.0+L2.*L3.*m2.*m3.*r2.*r3.*4.0+L1.*L3.*m2.*m3.*r2.*r3.*cos(th2).*4.0))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)-((tq1-L2.*g.*m3.*cos(th1+th2)-g.*m2.*r2.*cos(th1+th2)-L1.*g.*m2.*cos(th1)-L1.*g.*m3.*cos(th1)-g.*m1.*r1.*cos(th1)-g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth2.^2.*m3.*r3.*sin(th2+th3)+L1.*dth3.^2.*m3.*r3.*sin(th2+th3)+L1.*L2.*dth2.^2.*m3.*sin(th2)+L1.*dth2.^2.*m2.*r2.*sin(th2)+L2.*dth3.^2.*m3.*r3.*sin(th3)+L1.*dth1.*dth2.*m3.*r3.*sin(th2+th3).*2.0+L1.*dth1.*dth3.*m3.*r3.*sin(th2+th3).*2.0+L1.*dth2.*dth3.*m3.*r3.*sin(th2+th3).*2.0+L1.*L2.*dth1.*dth2.*m3.*sin(th2).*2.0+L1.*dth1.*dth2.*m2.*r2.*sin(th2).*2.0+L2.*dth1.*dth3.*m3.*r3.*sin(th3).*2.0+L2.*dth2.*dth3.*m3.*r3.*sin(th3).*2.0).*(Iz2.*Iz3-L2.^2.*L3.^2.*m3.^2-Iz3.*L2.^2.*m2-Iz2.*L3.^2.*m3+Iz3.*L2.^2.*m3+L2.^2.*L3.^2.*m2.*m3+L2.^2.*L3.*m3.^2.*r3.*2.0+Iz3.*L2.*m2.*r2.*2.0+Iz2.*L3.*m3.*r3.*2.0-L2.^2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.*L3.^2.*m3.^2.*cos(th2)+Iz3.*L1.*L2.*m3.*cos(th2)-L2.*L3.^2.*m2.*m3.*r2.*2.0-L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz3.*L1.*m2.*r2.*cos(th2)-L1.*L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3)+L1.*L2.*L3.*m3.^2.*r3.*cos(th2).*2.0-L1.*L3.^2.*m2.*m3.*r2.*cos(th2)+L2.*L3.*m2.*m3.*r2.*r3.*4.0+L1.*L3.*m2.*m3.*r2.*r3.*cos(th2).*2.0))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0);dth3;-((-tq3+g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth1.^2.*m3.*r3.*sin(th2+th3)+L2.*dth1.^2.*m3.*r3.*sin(th3)+L2.*dth2.^2.*m3.*r3.*sin(th3)+L2.*dth1.*dth2.*m3.*r3.*sin(th3).*2.0).*(Iz1.*Iz2+Iz1.*Iz3-L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m3.^2-L1.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*m1-Iz1.*L2.^2.*m2+Iz2.*L1.^2.*m2-Iz3.*L1.^2.*m1+Iz1.*L2.^2.*m3+Iz2.*L1.^2.*m3+Iz3.*L1.^2.*m2-Iz1.*L3.^2.*m3+Iz3.*L1.^2.*m3+L1.^2.*L2.^2.*m1.*m2-L1.^2.*L2.^2.*m1.*m3+L1.^2.*L3.^2.*m1.*m3-L1.^2.*L3.^2.*m2.*m3+L1.^2.*L2.*m2.^2.*r2.*2.0+L1.^2.*L3.*m3.^2.*r3.*2.0-L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*L1.*m1.*r1.*2.0+Iz3.*L1.*m1.*r1.*2.0+Iz1.*L2.*m2.*r2.*2.0+Iz1.*L3.*m3.*r3.*2.0-L1.^2.*m2.^2.*r2.^2.*cos(th2).^2+L1.^2.*L2.*m3.^2.*r3.*cos(th3).*2.0-L1.*L2.^2.*m1.*m2.*r1.*2.0+L1.*L2.^2.*m1.*m3.*r1.*2.0-L1.^2.*L2.*m1.*m2.*r2.*2.0-L1.*L3.^2.*m1.*m3.*r1.*2.0+L1.^2.*L2.*m2.*m3.*r2.*2.0-L1.^2.*L3.*m1.*m3.*r3.*2.0+L1.^2.*L3.*m2.*m3.*r3.*2.0+Iz1.*L2.*m3.*r3.*cos(th3).*2.0-L1.^2.*L2.*m3.^2.*r3.*cos(th2+th3).*cos(th2).*2.0-L1.^2.*L2.*m1.*m3.*r3.*cos(th3).*2.0+L1.^2.*L2.*m2.*m3.*r3.*cos(th3).*2.0+L1.*L2.*m1.*m2.*r1.*r2.*4.0+L1.*L3.*m1.*m3.*r1.*r3.*4.0-L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0-L1.^2.*m2.*m3.*r2.*r3.*cos(th2+th3).*cos(th2).*2.0+L1.*L2.*m1.*m3.*r1.*r3.*cos(th3).*4.0))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)+((-tq2+L2.*g.*m3.*cos(th1+th2)+g.*m2.*r2.*cos(th1+th2)+g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth1.^2.*m3.*r3.*sin(th2+th3)+L1.*L2.*dth1.^2.*m3.*sin(th2)+L1.*dth1.^2.*m2.*r2.*sin(th2)-L2.*dth3.^2.*m3.*r3.*sin(th3)-L2.*dth1.*dth3.*m3.*r3.*sin(th3).*2.0-L2.*dth2.*dth3.*m3.*r3.*sin(th3).*2.0).*(Iz1.*Iz3-L1.^2.*L3.^2.*m3.^2-Iz3.*L1.^2.*m1+Iz3.*L1.^2.*m2-Iz1.*L3.^2.*m3+Iz3.*L1.^2.*m3+L1.^2.*L3.^2.*m1.*m3-L1.^2.*L3.^2.*m2.*m3+L1.^2.*L3.*m3.^2.*r3.*2.0-L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2+Iz3.*L1.*m1.*r1.*2.0+Iz1.*L3.*m3.*r3.*2.0-L1.*L2.^2.*m3.^2.*r3.*cos(th2+th3)-L1.*L2.*L3.^2.*m3.^2.*cos(th2)+L1.^2.*L2.*m3.^2.*r3.*cos(th3)-Iz2.*L1.*m3.*r3.*cos(th2+th3)+Iz3.*L1.*L2.*m3.*cos(th2)-L1.*L3.^2.*m1.*m3.*r1.*2.0-L1.^2.*L3.*m1.*m3.*r3.*2.0+L1.^2.*L3.*m2.*m3.*r3.*2.0+Iz3.*L1.*m2.*r2.*cos(th2)+Iz1.*L2.*m3.*r3.*cos(th3)-L1.^2.*L2.*m3.^2.*r3.*cos(th2+th3).*cos(th2)-L1.*L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3)+L1.*L2.^2.*m3.^2.*r3.*cos(th2).*cos(th3)+L1.*L2.^2.*m2.*m3.*r3.*cos(th2+th3)+L1.*L2.*L3.*m3.^2.*r3.*cos(th2).*2.0-L1.*L3.^2.*m2.*m3.*r2.*cos(th2)-L1.^2.*L2.*m1.*m3.*r3.*cos(th3)+L1.^2.*L2.*m2.*m3.*r3.*cos(th3)+L1.*L3.*m1.*m3.*r1.*r3.*4.0-L1.^2.*m2.*m3.*r2.*r3.*cos(th2+th3).*cos(th2)-L1.*L2.*m2.*m3.*r2.*r3.*cos(th2+th3).*2.0+L1.*L2.*m1.*m3.*r1.*r3.*cos(th3).*2.0+L1.*L3.*m2.*m3.*r2.*r3.*cos(th2).*2.0+L1.*L2.*m2.*m3.*r2.*r3.*cos(th2).*cos(th3)))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)+(L1.*(-Iz2.*m3.*r3.*cos(th2+th3)+Iz3.*L2.*m3.*cos(th2)+Iz3.*m2.*r2.*cos(th2)-L2.^2.*m3.^2.*r3.*cos(th2+th3)-L2.*L3.^2.*m3.^2.*cos(th2)+L2.^2.*m2.*m3.*r3.*cos(th2+th3)+L2.*L3.*m3.^2.*r3.*cos(th2).*2.0-L3.^2.*m2.*m3.*r2.*cos(th2)-L2.*m3.^2.*r3.^2.*cos(th2+th3).*cos(th3)+L2.^2.*m3.^2.*r3.*cos(th2).*cos(th3)-L2.*m2.*m3.*r2.*r3.*cos(th2+th3).*2.0+L3.*m2.*m3.*r2.*r3.*cos(th2).*2.0+L2.*m2.*m3.*r2.*r3.*cos(th2).*cos(th3)).*(tq1-L2.*g.*m3.*cos(th1+th2)-g.*m2.*r2.*cos(th1+th2)-L1.*g.*m2.*cos(th1)-L1.*g.*m3.*cos(th1)-g.*m1.*r1.*cos(th1)-g.*m3.*r3.*cos(th1+th2+th3)+L1.*dth2.^2.*m3.*r3.*sin(th2+th3)+L1.*dth3.^2.*m3.*r3.*sin(th2+th3)+L1.*L2.*dth2.^2.*m3.*sin(th2)+L1.*dth2.^2.*m2.*r2.*sin(th2)+L2.*dth3.^2.*m3.*r3.*sin(th3)+L1.*dth1.*dth2.*m3.*r3.*sin(th2+th3).*2.0+L1.*dth1.*dth3.*m3.*r3.*sin(th2+th3).*2.0+L1.*dth2.*dth3.*m3.*r3.*sin(th2+th3).*2.0+L1.*L2.*dth1.*dth2.*m3.*sin(th2).*2.0+L1.*dth1.*dth2.*m2.*r2.*sin(th2).*2.0+L2.*dth1.*dth3.*m3.*r3.*sin(th3).*2.0+L2.*dth2.*dth3.*m3.*r3.*sin(th3).*2.0))./(Iz1.*Iz2.*Iz3-Iz2.*Iz3.*L1.^2.*m1-Iz1.*Iz3.*L2.^2.*m2+Iz2.*Iz3.*L1.^2.*m2-Iz1.*Iz2.*L3.^2.*m3+Iz1.*Iz3.*L2.^2.*m3+Iz2.*Iz3.*L1.^2.*m3-L1.^2.*L2.^2.*L3.^2.*m3.^3-Iz3.*L1.^2.*L2.^2.*m2.^2-Iz1.*L2.^2.*L3.^2.*m3.^2-Iz2.*L1.^2.*L3.^2.*m3.^2+Iz3.*L1.^2.*L2.^2.*m3.^2-Iz3.*L1.^2.*L2.^2.*m3.^2.*cos(th2).^2+Iz2.*Iz3.*L1.*m1.*r1.*2.0+Iz1.*Iz3.*L2.*m2.*r2.*2.0+Iz1.*Iz2.*L3.*m3.*r3.*2.0-Iz3.*L1.^2.*m2.^2.*r2.^2.*cos(th2).^2-Iz1.*L2.^2.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.*m3.^3.*r3.*2.0-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).^2+L1.^2.*L2.^2.*L3.^2.*m3.^3.*cos(th2).^2-L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th3).^2+L1.^2.*L2.^2.*L3.^2.*m1.*m3.^2+L1.^2.*L2.^2.*L3.^2.*m2.^2.*m3+Iz3.*L1.^2.*L2.^2.*m1.*m2+Iz2.*L1.^2.*L3.^2.*m1.*m3-Iz3.*L1.^2.*L2.^2.*m1.*m3+Iz1.*L2.^2.*L3.^2.*m2.*m3-Iz2.*L1.^2.*L3.^2.*m2.*m3+Iz3.*L1.^2.*L2.*m2.^2.*r2.*2.0+Iz1.*L2.^2.*L3.*m3.^2.*r3.*2.0+Iz2.*L1.^2.*L3.*m3.^2.*r3.*2.0-Iz2.*L1.^2.*m3.^2.*r3.^2.*cos(th2+th3).^2-L1.^2.*L2.^2.*L3.*m3.^3.*r3.*cos(th2).^2.*2.0-L1.^2.*L2.^2.*L3.^2.*m1.*m2.*m3-L1.*L2.^2.*L3.^2.*m1.*m3.^2.*r1.*2.0-L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*2.0-L1.^2.*L2.*L3.^2.*m2.^2.*m3.*r2.*2.0-L1.^2.*L2.^2.*L3.*m1.*m3.^2.*r3.*2.0-L1.^2.*L2.^2.*L3.*m2.^2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th2+th3).^2-Iz3.*L1.*L2.^2.*m1.*m2.*r1.*2.0-Iz2.*L1.*L3.^2.*m1.*m3.*r1.*2.0+Iz3.*L1.*L2.^2.*m1.*m3.*r1.*2.0-Iz3.*L1.^2.*L2.*m1.*m2.*r2.*2.0-Iz1.*L2.*L3.^2.*m2.*m3.*r2.*2.0-Iz2.*L1.^2.*L3.*m1.*m3.*r3.*2.0+Iz3.*L1.^2.*L2.*m2.*m3.*r2.*2.0-Iz1.*L2.^2.*L3.*m2.*m3.*r3.*2.0+Iz2.*L1.^2.*L3.*m2.*m3.*r3.*2.0+L1.^2.*L2.^2.*m1.*m3.^2.*r3.^2.*cos(th3).^2+L1.^2.*L3.^2.*m2.^2.*m3.*r2.^2.*cos(th2).^2-L1.^2.*L2.^2.*m2.*m3.^2.*r3.^2.*cos(th3).^2-L1.*L2.^2.*m1.*m3.^2.*r1.*r3.^2.*cos(th3).^2.*2.0-L1.^2.*L3.*m2.^2.*m3.*r2.^2.*r3.*cos(th2).^2.*2.0-Iz3.*L1.^2.*L2.*m2.*m3.*r2.*cos(th2).^2.*2.0+L1.^2.*L2.^2.*m3.^3.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.^2.*L3.^2.*m1.*m2.*m3.*r1.*2.0+L1.^2.*L2.*L3.^2.*m1.*m2.*m3.*r2.*2.0+L1.^2.*L2.^2.*L3.*m1.*m2.*m3.*r3.*2.0+L1.*L2.^2.*L3.*m1.*m3.^2.*r1.*r3.*4.0+L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*4.0+L1.^2.*L2.*L3.*m2.^2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).^2.*2.0+L1.^2.*L2.*L3.^2.*m2.*m3.^2.*r2.*cos(th2).^2.*2.0+Iz3.*L1.*L2.*m1.*m2.*r1.*r2.*4.0+Iz2.*L1.*L3.*m1.*m3.*r1.*r3.*4.0+Iz1.*L2.*L3.*m2.*m3.*r2.*r3.*4.0-L1.*L2.*L3.^2.*m1.*m2.*m3.*r1.*r2.*4.0-L1.*L2.^2.*L3.*m1.*m2.*m3.*r1.*r3.*4.0-L1.^2.*L2.*L3.*m1.*m2.*m3.*r2.*r3.*4.0-L1.^2.*L2.*L3.*m2.*m3.^2.*r2.*r3.*cos(th2).^2.*4.0+L1.^2.*L2.*m2.*m3.^2.*r2.*r3.^2.*cos(th2+th3).*cos(th2).*cos(th3).*2.0+L1.*L2.*L3.*m1.*m2.*m3.*r1.*r2.*r3.*8.0)];
function Jac = get_Jacobian(th1, th2, th3)

global  L1 L2 L3;

J11 = -L1*sin(th1) - L2*sin(th1+th2) - L3*sin(th1+th2+th3);
J12 = -L2*sin(th1+th2) - L3*sin(th1+th2+th3);
J13 = -L3*sin(th1+th2+th3);

J21 = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
J22 = L2*cos(th1+th2) + L3*cos(th1+th2+th3);
J23 = L3*cos(th1+th2+th3);

Jac = [J11, J12, J13;
       J21, J22, J23];

end
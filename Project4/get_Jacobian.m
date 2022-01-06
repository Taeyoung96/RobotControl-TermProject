function Jac = get_Jacobian(th1, th2)

global  L1 L2;

Jac = [-L1*sin(th1)- L2*sin(th1+th2),  -L2*sin(th1+th2);
        L1*cos(th1)+ L2*cos(th1+th2),   L2*cos(th1+th2)];

end


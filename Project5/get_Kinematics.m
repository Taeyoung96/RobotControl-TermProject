function [X] = get_Kinematics(th1,th2, th3)
    global L1 L2 L3
    
    
    X = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    Y = L1*sin(th1) + L2*sin(th1+th2) + L3*sin(th1+th2+th3);
    
    X = [X;Y];
end
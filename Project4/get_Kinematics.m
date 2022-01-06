function [X] = get_Kinematics(th1,th2)
    global L1 L2
    
    
    X = L1*cos(th1) + L2*cos(th1+th2);
    Y = L1*sin(th1) + L2*sin(th1+th2);
    
    X = [X;Y];
end


function dxdt = one_link(t,x)

global m L I tq g

r = (1/3)*L;



                                         
dxdt = [ x(2);                        % dx
        (tq - m*r*g*sin(x(1)))/(I)];  % ddx
    
    
    
    %(tq - m*g*r*cos(x(1))) / (- m*L^2 + 2*m*r*L+I)
    
    %(tq - m*r*g*sign(x(1)))/(I)
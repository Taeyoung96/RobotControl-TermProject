function T = HT(th, d, a, aI)

T = [ cos(th), -cos(aI)*sin(th),   sin(aI)*sin(th),    a*cos(th);
      sin(th),  cos(aI)*cos(th),  -sin(aI)*cos(th),    a*sin(th);
      0      ,  sin(aI)        ,   cos(aI)        ,    d        ;
      0      ,  0              ,   0              ,    1       ]; 


end




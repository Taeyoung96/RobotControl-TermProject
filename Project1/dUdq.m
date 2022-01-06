function DQ=dUdq(i,j,k)

global T01 T02 T03 T12 T23 T13
global U11 U12 U13 U21 U22 U23 U31 U32 U33
global J1 J2 J3
global Q1 Q2 Q3

T00 = eye(4);
T11 = eye(4);
T22 = eye(4);
T33 = eye(4);

if (i>=k)&(k>=j)
    cmd = sprintf("DQ = T0%d * Q%d * T%d%d * Q%d * T%d%d;",j-1,j,j-1,k-1,k,k-1,i);
elseif (i>=j)&(j>k)
    cmd = sprintf("DQ = T0%d * Q%d * T%d%d * Q%d * T%d%d;",k-1,k,k-1,j-1,j,j-1,i);
elseif (i<j)|(i<k)
    cmd = sprintf("DQ = zeros(4,4);");
end
eval(cmd)
end
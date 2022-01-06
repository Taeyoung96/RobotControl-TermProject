function dxdt = TwoLinkRobot(t,x)
 
global I1; global I2; global L1; global L2; global Im1; global Im2; global m1; global m2;
global g; global r1; global r2; global Fs1; global Fs2; global Fv1; global Fv2; global tau1; global tau2;
 
q1 = x(1); q1dot = x(2);
q2 = x(3); q2dot = x(4);
 
M_11 = I1 + I2 + m2*L1*L1 + 2*m2*r2*L1*cos(q2) + Im1;
M_12 = I2 + m2*r2*L1*cos(q2);
M_21 = I2 + m2*r2*L1*cos(q2);
M_22 = I2 + Im2;
C_11 = -m2*r2*L1*sin(q2)*q2dot;
C_12 = -m2*r2*L1*sin(q2)*(q1dot+q2dot);
C_21 = m2*r2*L1*sin(q2)*q1dot;
C_22 = 0;
g_1 = -m1*r1*g*cos(q1) - m2*L1*g*cos(q1) - m2*r2*g*cos(q1+q2);
g_2 = -m2*r2*g*cos(q1+q2);
d_1 = Fs1*sign(q1dot) + Fv1*q1dot;
d_2 = Fs2*sign(q2dot) + Fv2*q2dot;
 
C = [C_11 C_12 ; C_21 C_22];
G = [g_1 ; g_2];
D = [d_1 ; d_2];
M = [M_11, M_12 ; M_21, M_22];
R = C * [q1dot ; q2dot] + G + D;
 
tau = [tau1 ; tau2] - R;
 
vdot = M\tau;
dxdt = [q1dot ; vdot(1) ; q2dot ; vdot(2)];
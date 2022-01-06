%% Using the Lagrangian function, derive 3 DOF robot and perform a free fall simulation. Discuss the results.
% Transoformation matrix

clc
clear all
close all

syms L1 L2 L3 m1 m2 m3 Ic1 Ic2 Ic3 
syms Im1 Im2 Im3 r1 r2 r3 Iz1 Iz2 Iz3
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3

I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;

I2xx = 0;
I2yy = Iz2;
I2zz = Iz2;

I3xx = 0;
I3yy = Iz3;
I3zz = Iz3;

d1 = 0;     d2 = 0;     d3 = 0;
a1 = L1;    a2 = L2;    a3 = L3;
aI1 = 0;    aI2 = 0;    aI3 = 0;

DH = [th1, d1, a1, aI1;
      th2, d2, a2, aI2;
      th3, d3, a3, aI3];
  
global T01 T02 T03 T12 T23 T13

T01 = HT(th1, d1, a1, aI1);
T12 = HT(th2, d2, a2, aI2);
T23 = HT(th3, d3, a3, aI3);

T02 = T01*T12;
T03 = T01*T12*T23;
T13 = T12*T23;

%% U matrix

Qr = [0 -1  0  0 ;
      1  0  0  0 ;
      0  0  0  0 ;
      0  0  0  0];

global Q1 Q2 Q3
  
Q1 = Qr;
Q2 = Qr;
Q3 = Qr;

global U11 U12 U13 U21 U22 U23 U31 U32 U33

U11 = Q1*T01;       % i=1 j=1
U12 = zeros(4,4);   % i=1 j=2
U13 = zeros(4,4);   % i=1 j=3

U21 = Q1*T02;       % i=2 j=1
U22 = T01*Q2*T12;   % i=2 j=2
U23 = zeros(4,4);   % i=2 j=3

U31 = Q3*T03;       % i=3 j=1
U32 = T01*Q2*T13;   % i=3 j=2
U33 = T02*Q3*T23;   % i=3 j=3

%% Pseudo-inverse
global J1 J2 J3

J1 = [1/2 * (-I1xx + I1yy + I1zz) 0                           0                           -m1*(L1-r1);
      0                           1/2 * (I1xx - I1yy + I1zz)  0                           0;
      0                           0                           1/2 * (I1xx + I1yy - I1zz)  0;
      -m1*(L1-r1)                 0                           0                           m1];

J2 = [1/2 * (-I2xx + I2yy + I2zz) 0                           0                           -m2*(L2-r2);
      0                           1/2 * (I2xx - I2yy + I2zz)  0                           0;
      0                           0                           1/2 * (I2xx + I2yy - I2zz)  0;
      -m2*(L2-r2)                 0                           0                           m2];

J3 = [1/2 * (-I3xx + I3yy + I3zz) 0                           0                           -m3*(L3-r3);
      0                           1/2 * (I3xx - I3yy + I3zz)  0                           0;
      0                           0                           1/2 * (I3xx + I3yy - I3zz)  0;
      -m3*(L3-r3)                 0                           0                           m3];
  
%% Inertia matrix

n = 3;

for i=1:n
    for k=1:n
        M(i,k) = Inertia(i,k,n);
    end
end

%% dUdq

n = 3;

for i=1:n
    for j=1:n
        for k=1:n
            cmd = sprintf('U%d%d%d = dUdq(i,j,k);',i,j,k);
            eval(cmd);
        end
    end
end

%% Lagrangian function 
                                    %                 dth(k)*dth(m)
h111 = trace(U111*J1*U11.')...      % i=1 k=1 m=1 j=1   
      +trace(U211*J2*U21.')...      % i=1 k=1 m=1 j=2   
      +trace(U311*J3*U31.');        % i=1 k=1 m=1 j=3   dth1*dth1
  
h112 = trace(U212*J2*U21.')...      % i=1 k=1 m=2 j=2   
      +trace(U312*J3*U31.');        % i=1 k=1 m=2 j=3   dth1*dth2

h113 = trace(U313*J3*U31.');        % i=1 k=1 m=3 j=3   dth1*dth3
  

h121 = trace(U221*J2*U21.')...      % i=1 k=2 m=1 j=2   
      +trace(U321*J3*U31.');        % i=1 k=2 m=1 j=3   dth2*dth1

h122 = trace(U222*J2*U21.')...      % i=1 k=2 m=2 j=2   
      +trace(U322*J3*U31.');        % i=1 k=2 m=2 j=3   dth2*dth2
  
h123 = trace(U323*J3*U31.');        % i=1 k=2 m=3 j=3   dth2*dth3

h131 = trace(U331*J3*U31.');        % i=1 k=3 m=1 j=3   dth3*dth1

h132 = trace(U332*J3*U31.');        % i=1 k=3 m=2 j=3   dth3*dth2

h133 = trace(U333*J3*U31.');        % i=1 k=3 m=3 j=3   dth3*dth3


h1 = (dth1*dth1)*(h111) + (dth1*dth2)*(h112) + (dth1*dth3)*(h113)...
    +(dth2*dth1)*(h121) + (dth2*dth2)*(h122) + (dth2*dth3)*(h123)...
    +(dth3*dth1)*(h131) + (dth3*dth2)*(h132) + (dth3*dth3)*(h133);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    %                 dth(k)*dth(m)
h211 = trace(U211*J2*U22.')...      % i=2 k=1 m=1 j=2   
      +trace(U311*J3*U32.');        % i=2 k=1 m=1 j=3   dth1*dth1

h212 = trace(U212*J2*U22.')...      % i=2 k=1 m=2 j=2
      +trace(U312*J3*U32.');        % i=2 k=1 m=2 j=3   dth1*dth2

h213 = trace(U313*J3*U32.');        % i=2 k=1 m=3 j=3   dth1*dth3
  
h221 = trace(U221*J2*U22.')...      % i=2 k=2 m=1 j=2   dth2*dth1
      +trace(U321*J3*U32.');        % i=2 k=2 m=1 j=3

h222 = trace(U222*J2*U22.')...      % i=2 k=2 m=2 j=2   dth2*dth2
      +trace(U322*J3*U32.');        % i=2 k=2 m=2 j=3

h223 = trace(U323*J3*U32.');        % i=2 k=2 m=3 j=3   dth2*dth3

h231 = trace(U331*J3*U32.');        % i=2 k=3 m=1 j=3   dth3*dth1

h232 = trace(U332*J3*U32.');        % i=2 k=3 m=2 j=3   dth3*dth2

h233 = trace(U333*J3*U32.');        % i=2 k=3 m=3 j=3   dth3*dth3

h2 = (dth1*dth1)*(h211) + (dth1*dth2)*(h212) + (dth1*dth3)*(h213)...
    +(dth2*dth1)*(h221) + (dth2*dth2)*(h222) + (dth2*dth3)*(h223)...
    +(dth3*dth1)*(h231) + (dth3*dth2)*(h232) + (dth3*dth3)*(h233);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    %                 dth(k)*dth(m)
h311 = trace(U311*J3*U33.');        % i=3 k=1 m=1 j=3   dth1*dth1

h312 = trace(U312*J3*U33.');        % i=3 k=1 m=2 j=3   dth1*dth2

h313 = trace(U313*J3*U33.');        % i=3 k=1 m=3 j=3   dth1*dth3

h321 = trace(U321*J3*U33.');        % i=3 k=2 m=1 j=3   dth2*dth1

h322 = trace(U322*J3*U33.');        % i=3 k=2 m=2 j=3   dth2*dth2

h323 = trace(U323*J3*U33.');        % i=3 k=2 m=3 j=3   dth2*dth3

h331 = trace(U331*J3*U33.');        % i=3 k=3 m=1 j=3   dth3*dth1

h332 = trace(U332*J3*U33.');        % i=3 k=3 m=2 j=3   dth3*dth2

h333 = trace(U333*J3*U33.');        % i=3 k=3 m=3 j=3   dth3*dth3

h3 = (dth1*dth1)*(h311) + (dth1*dth2)*(h312) + (dth1*dth3)*(h313)...
    +(dth2*dth1)*(h321) + (dth2*dth2)*(h322) + (dth2*dth3)*(h323)...
    +(dth3*dth1)*(h331) + (dth3*dth2)*(h332) + (dth3*dth3)*(h333);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H matrix

h = simplify([h1; h2; h3]);

%% Gravity term 도출

syms g
r11 = [-(L1-r1); 0; 0; 1];  %무게중심까지의 거리
r22 = [-(L2-r2); 0; 0; 1];
r33 = [-(L3-r3); 0; 0; 1];

gv = [0 -g 0 0];

G1 = -( m1*gv*U11*r11...   % i=1 j=1
       +m2*gv*U21*r22...   % i=1 j=2
       +m3*gv*U31*r33);    % i=1 j=3
   
   
G2 = -( m2*gv*U22*r22...   % i=2 j=2
       +m3*gv*U32*r33);    % i=2 j=3
   
G3 = -(m3*gv*U33*r33);     % i=3 j=3

G = simplify([G1;G2;G3]);

%% 강의자료 p14 - 확장

syms tau1 tau2 tau3

DDTH = inv(M)*([tau1; tau2; tau3] -h -G);  % theta 1 dot과 theta 2dot에 의해 만들어진 행렬


%% Y와 dydt 구하기

dydt = simplify([dth1; DDTH(1); dth2; DDTH(2); dth3; DDTH(3)]);
%matlabFunction(dydt,'file','three_links.m','Optimize',false); % 위 dydt를 함수로 만들기

%% 그래프 띄우기
clear all
close all

% Initialize
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tau1 tau2 tau3

L1 = 0.5;   L2 = 0.5;   L3 = 0.5;
r1 = 0.1;   r2 = 0.1;   r3 = 0.1;
m1 = 0.2;   m2 = 0.2;   m3 = 0.2;

Iz1 = 0.05; Iz2 = 0.05; Iz3 = 0.05;

g = 9.806;

dt = 0.02; ft = 5;

q1 = -pi/2; dq1 = 0;
q2 = pi/4;  dq2 = 0;
q3 = pi/4;  dq3 = 0;

data = [];

n = 1; 

FG = figure('Position',[300 300 600 600],'Color',[1 1 1])
AX = axes('parent',FG);

hold on
grid on
axis([-2.0 2.0 -2.0 2.0]);

x1 = L1*cos(q1);
y1 = L1*sin(q1);
Px1 = [0 x1];
Py1 = [0 y1];

x2 = L2*cos(q1+q2);
y2 = L2*sin(q1+q2); 
Px2 = [x1 x1+x2];
Py2 = [y1 y1+y2];

x3 = L3*cos(q1+q2+q3);
y3 = L3*sin(q1+q2+q3);
Px3 = [x1+x2 x1+x2+x3];
Py3 = [y1+y2 y1+y2+y3];

p1 = plot(Px1,Py1,'-ob','Linewidth',3);
p2 = plot(Px2,Py2,'-or','Linewidth',3);
p3 = plot(Px3,Py3,'-og','Linewidth',3);


for cnt=0:dt:ft
    
    tau1 = 0.0;
    tau2 = 0.0;
    tau3 = 0.0;
    
    [t,y] = ode45('three_links',[0 dt],[q1; dq1; q2; dq2; q3; dq3]);
    
    index = length(y);
    
    q1  = y(index,1);
    dq1 = y(index,2);
    q2  = y(index,3);
    dq2 = y(index,4);
    q3  = y(index,5);
    dq3 = y(index,6);
    
    x1 = L1*cos(q1);
    y1 = L1*sin(q1);
    Px1 = [0 x1];
    Py1 = [0 y1];
    
    x2 = L2*cos(q1+q2);
    y2 = L2*sin(q1+q2); 
    Px2 = [x1 x1+x2];
    Py2 = [y1 y1+y2];
    
    x3 = L3*cos(q1+q2+q3);
    y3 = L3*sin(q1+q2+q3);
    Px3 = [x1+x2 x1+x2+x3];
    Py3 = [y1+y2 y1+y2+y3];
    
    n = n+1;
    
    cmd = sprintf('Time : %2.2f',cnt);
    clc
    disp(cmd)
    
    if rem(n,2)==0
        set(p1,'XData',Px1,'YData',Py1)
        set(p2,'XData',Px2,'YData',Py2)
        set(p3,'XData',Px3,'YData',Py3)
        drawnow
    end
end
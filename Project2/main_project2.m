%% Case 1: Kalman filter based error estimation

% Initialization
clc
clear all
close all

global I1; 
global I2; 
global L1; 
global L2; 
global Im1; 
global Im2; 
global m1; 
global m2;
global g; 
global r1; 
global r2; 
global Fs1; 
global Fs2; 
global Fv1; 
global Fv2; 
global tau1; 
global tau2;

I1 = 0.05; 
I2 = 0.05; 
Im1 = 0.05; 
Im2 = 0.05; 
m1 = 0.2; 
m2 = 0.2; 
g = 9.806; 
Fs1 = 0.1; 
Fs2 = 0.1; 
Fv1 = 0.1; 
Fv2 = 0.1; 
L1 = 0.5; 
L2 = 0.5; 
r1 = 0.1; 
r2 = 0.1; 

dT = 0.002;
ft = 5;
q1 = 0;
q1dot = 0;
q2 = 0;
q2dot = 0; 
n = 1;

w1_int = zeros(2,10);
u = 0;           
P = eye(10,10);  
R = eye(2,2);

theta = [0;0;0;0;0;0;0;0;0;0];
Yn = zeros(10);
un = 0;

%% Estimate unknwon parameters with regressor

for t = 0 : dT : ft
   clc
   disp(t)
    
   %%%% Dynamics
   tau1 = (sin(t) + cos(10*t));  % 임의의 연속적 입력토크 tau 생성
   tau2 = (sin(10*t) + cos(10*t));
   tau = [tau1 ; tau2];
   [st,x] = ode45('TwoLinkRobot',[0,dT],[q1; q1dot; q2; q2dot]); % Dynmaics 실행
   index = size(x); 
   
   q1    =  x(index(1), 1) ; % 결과 q값 저장
   q1dot =  x(index(1), 2);  % 결과 dq값 저장
   q2    =  x(index(1), 3);
   q2dot =  x(index(1), 4); 
   
   %%%% Regressor
   w1_int = w1_int + [ 0  0   0                               g*cos(q1)   g*cos(q1+q2)    0   -sign(q1dot)    0               -q1dot  0     ; 
                       0  0   -sin(q2)*(q1dot+q2dot)*q1dot    0           g*cos(q1+q2)    0   0               -sign(q2dot)    0       -q2dot] * dT;
                   
                   
   w2     =          [q1dot    q1dot+q2dot     2*cos(q2)*(q1dot+q2dot)     0   0   0       0   0   0   0 ; 
                        0      q1dot+q2dot     cos(q2)*q1dot               0   0   q2dot   0   0   0   0 ];
   Y = w2 - w1_int;
   u = u + tau*dT;
   
   %%%% Case 1 : Kalman filter based parameter estimation algorithm
    P = P - P * Y' * inv(R + Y * P * Y') * Y * P;
    K = P * Y';
    theta = theta + K*(u - Y*theta);

   %%%% Data save   
   save_time(n,:) = t;       % 시간 저장
   save_theta(n,:) = theta;  % 결과 저장
   n = n + 1;    
   
end
%% Result Graph Visualize

disp(theta) % theta 최종 결과  command 창에서 확인

figure(1)
hold on
grid on
plot(save_time, save_theta)
xlabel('Time[s]');
ylabel('Theta(Unknown Parameter)');
legend('I1+m2*L1*L1+Im1','I2','m2*r2*L1','m1*r1+m2*L1','m2*r2','Im2','Fs1','Fs2','Fv1','Fv2');
xlim ([0 5])

%% Case 2 : Error minimization algoritm

% Initialization
clc
clear all
close all

global I1; 
global I2; 
global L1; 
global L2; 
global Im1; 
global Im2; 
global m1; 
global m2;
global g; 
global r1; 
global r2; 
global Fs1; 
global Fs2; 
global Fv1; 
global Fv2; 
global tau1; 
global tau2;

I1 = 0.05; 
I2 = 0.05; 
Im1 = 0.05; 
Im2 = 0.05; 
m1 = 0.2; 
m2 = 0.2; 
g = 9.806; 
Fs1 = 0.1; 
Fs2 = 0.1; 
Fv1 = 0.1; 
Fv2 = 0.1; 
L1 = 0.5; 
L2 = 0.5; 
r1 = 0.1; 
r2 = 0.1; 

dT = 0.002;
ft = 5;
q1 = 0;
q1dot = 0;
q2 = 0;
q2dot = 0; 
n = 1;

w1_int = zeros(2,10);
u = 0;           
P = eye(10,10);  
R = eye(2,2);

theta = [0;0;0;0;0;0;0;0;0;0];
Yn = zeros(10);
un = 0;

%% Estimate unknwon parameters with regressor
for t = 0 : dT : ft
   clc
   disp(t)
    
   %%%% Dynamics
   tau1 = (sin(t) + cos(10*t));  % 임의의 연속적 입력토크 tau 생성
   tau2 = (sin(10*t) + cos(10*t));
   tau = [tau1 ; tau2];
   [st,x] = ode45('TwoLinkRobot',[0,dT],[q1; q1dot; q2; q2dot]); % Dynmaics 실행
   index = size(x); 
   
   q1    =  x(index(1), 1) ; % 결과 q값 저장
   q1dot =  x(index(1), 2);  % 결과 dq값 저장
   q2    =  x(index(1), 3);
   q2dot =  x(index(1), 4); 
   
   %%%% Regressor
   w1_int = w1_int + [ 0  0   0                               g*cos(q1)   g*cos(q1+q2)    0   -sign(q1dot)    0               -q1dot  0     ; 
                       0  0   -sin(q2)*(q1dot+q2dot)*q1dot    0           g*cos(q1+q2)    0   0               -sign(q2dot)    0       -q2dot] * dT;
                   
                   
   w2     =          [q1dot    q1dot+q2dot     2*cos(q2)*(q1dot+q2dot)     0   0   0       0   0   0   0 ; 
                        0      q1dot+q2dot     cos(q2)*q1dot               0   0   q2dot   0   0   0   0 ];
   Y = w2 - w1_int;
   u = u + tau*dT;
   
  
   %%%% Case 2 : Error minimization algoritm
  Yn = Yn + Y' * Y;
  un = un + Y' * u;
  theta = inv(Yn)*un;
   
   %%%% Data save   
   save_time(n,:) = t;       % 시간 저장
   save_theta(n,:) = theta;  % 결과 저장
   n = n + 1;    
   
end
%%
%%%% 결과 그래프 도출
disp(theta) % theta 최종 결과  command 창에서 확인

figure(1)
hold on
grid on
plot(save_time, save_theta)
xlabel('Time[s]');
ylabel('Theta(Unknown Parameter)');
legend('I1+m2*L1*L1+Im1','I2','m2*r2*L1','m1*r1+m2*L1','m2*r2','Im2','Fs1','Fs2','Fv1','Fv2');
xlim ([0 5])
%Initialization
clc,clear all 
close all

%% Set Simulation Parameters      EXP.3 ->   2. Joint Space PD CTM Controller
    %Draw flag
    flag_Simul       =1;
    flag_Draw        =1;
    flag_Draw_Robot = 1;
    flag_Draw_Graph = 1;

    %Global Variable
        global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tq1 tq2 tq3;

    %simulation Parameters
        delta_t     = 0.005;    %[sec]      :Sampling Time
        start_t     = 0.000;    %[sec]      :Start Time
        finish_t    = 6.000;    %[sec]      : End Time

        g           =9.8148;    %[m/s^2]    : Gravitational Acceleration

        %Robot Parameters
        m1  = 0.2;   m2  = 0.2;    m3 = 0.2;              %[kg]       :Link Mass
        L1  = 0.50;  L2  = 0.50;   L3 = 0.5;              %[m]       :Link Length
        r1  = 0.10;  r2  = 0.10;   r3 = 0.1;              %[m]       :Link Centor of mass
        Iz1 = 0.05;  Iz2 = 0.05;   Iz3 = 0.05;            %[kgm^2]       :Link Inertia

        init_q1 = -pi/2; init_q2 = 0; init_q3 = 0;       
        init_dq1 = 0.0;  init_dq2 = 0.0; init_dq3 = 0.0;       
        init_ddq1 = 0.0;  init_ddq2 = 0.0; init_ddq3 = 0.0;       

        init_q           = [init_q1; init_q2; init_q3];           %[rad]      : Init Joint Angle 
        init_dq          = [init_dq1; init_dq2; init_dq3];        %[rad/s]    : Init Angular Velocity
        init_ddq          = [init_ddq1; init_ddq2; init_ddq3];    %[rad/s^2]    : Init Angular Acceleration
        
        q = init_q;                         %[rad]      : Current Joint Angle
        dq = init_dq;                       %[rad/s]    : Current Angular Velocity
        ddq = init_ddq;                     %[rad/s]    : Current Angular Velocity
        
        q_d_int=0;
        
        % Target Position Parameters
        q_d         =   init_q;             %[rad]          : Target Joint Angle
        dq_d        =   init_dq;            %[rad/s]        : Target Angular Velocity
        ddq_d       =   init_ddq;                  %[rad/s^2]      : Target Angular Acceleration
        
        tq1 = 0.0;              % [Nm]  : Control Torque1
        tq2 = 0.0;              % [Nm]  : Control Torque2
        tq3 = 0.0;              % [Nm]  : Control Torque3
        tq = [tq1; tq2; tq3];              % [Nm]  : Control Torque
             
        %Controller Gain
        zeta        =   1;    %critical damped system
        Wn          =   20;         %[rad/s]        : natural frequency
        Kp          =   Wn^2;       %[Nm/rad]       : propotional Gain
        Kv          =   2*zeta*Wn;  %[Nm*s/rad]:    : Derivative Gain
        
        Ki          = 1000;
            
%% Simulation
if(flag_Simul == 1)
    %Simulation
        n = 1;
        sin_t = 0;
        for(time = start_t:delta_t:finish_t)
            %Set Target Position   
            if(time < 1.0)
                q_d         = init_q;
                dq_d        = init_dq;
                ddq_d       = init_ddq;
            else
                if(q_d(1) < 0*pi/180) % 현재 Desired 값이 90도 보다 작으면
                    q_d(1) = q_d(1) + (45*(pi/180)/2)*delta_t;      
                else
                    q_d(1) = 0*pi/180;
                end
                        dq_d(1)    = 30*pi/180;    %(q_d - simul_q_d(n-1))/delta_t;          %simul_q_d = pre_q_d
                        ddq_d(1)   = (dq_d(1) - simul_dq_d1(n-1))/delta_t;        %simul_dq_d = pre_dq_d
            end
            
            % Error Sum
            q_d_int = q_d_int + (q_d-q)*delta_t;
                    
            %Get Dynamics
            I = get_Inertia(q(2),q(3));                  %inertia
            C = get_Coriollis(q(2),q(3), dq(1),dq(2),dq(3));         %coriollis
            G = get_Gravity(q(1),q(2),q(3));             % gravity

            % Controller
            %u        = ddq_d + Kp*(q_d - q);                             % P Controller
            %u        = ddq_d + Kv*(dq_d - dq) + Kp*(q_d - q);            % PD Controller
            u        = ddq_d + Kv*(dq_d - dq) + Kp*(q_d - q)+Ki*q_d_int; %PID Controller
            tq_ctrl = I*u + C*dq + G*1;         

             %Robot Model
                % Inverse Dynamics
                tq = tq_ctrl;
                tq1 = tq(1);    tq2 = tq(2);     tq3 = tq(3);
              
                [t,y] = ode45('three_links', [0 delta_t], [q(1); dq(1); q(2); dq(2); q(3); dq(3)]);
                index = length(y);
                q =  [y(index, 1); y(index,3); y(index,5)];
                dq = [y(index,2);  y(index,4); y(index,6)];

                %save Data
                simul_time(n)   = time;          %[sec]
                
                % Current
                simul_q1(n)      = q(1);         %[rad]
                simul_q2(n)      = q(2);         %[rad]
                simul_q3(n)      = q(3);         %[rad]
                
                simul_dq1(n)     =  dq(1);       % [rad/s]
                simul_dq2(n)      = dq(2);       % [rad/s]
                simul_dq3(n)      = dq(3);       % [rad/s]
                
                % Desired
                simul_q_d1(n)     = q_d(1);      %[rad]
                simul_q_d2(n)     = q_d(2);      %[rad]
                simul_q_d3(n)     = q_d(3);      %[rad]
                
                simul_dq_d1(n)     =  dq_d(1);       % [rad/s]
                simul_dq_d2(n)      = dq_d(2);       % [rad/s]
                simul_dq_d3(n)      = dq_d(3);       % [rad/s]
              
                n                 = n+ 1;
        end
end
            
    
%% Simulation Result Graph
if(flag_Draw == 1)
    fontsize           =20;
    font_size_title    =25;
    font_size_label    =20;
    linewidth_current   =3;
    linewidth_target    =5;

    if(flag_Draw_Robot == 1)
        % Draw Robot
            x1 = L1*cos(init_q1);               % [m]  : Joint 1 X-axis position
            y1 = L1*sin(init_q1);               % [m]  : Joint 1 Y-axis position
            x2 = L2*cos(init_q1+init_q2);       % [m]  : Joint 2 X-axis position
            y2 = L2*sin(init_q1+init_q2);       % [m]  : Joint 2 Y-axis position
            x3 = L3*cos(init_q1+init_q2+init_q3);       % [m]  : Joint 2 X-axis position
            y3 = L3*sin(init_q1+init_q2+init_q3);       % [m]  : Joint 2 Y-axis position

            FG1 = figure('Position',[200 300 700 700], 'Color', [1 1 1]);
                AX = axes('parent',FG1); hold on

                Px1 = [0  x1];
                Py1 = [0  y1];
                Px2 = [x1 x1+x2];
                Py2 = [y1 y1+y2];
                Px3 = [x1+x2 x1+x2+x3];
                Py3 = [y1+y2 y1+y2+y3];

                p1 = plot(Px1,Py1, '-ob', 'Linewidth', linewidth_current);
                p2 = plot(Px2,Py2, '-or', 'Linewidth', linewidth_current);
                p3 = plot(Px3,Py3, '-og', 'Linewidth', linewidth_current);

                axis([-2.0 2.0 -2.0 2.0]);
                grid on
            xlabel('X-axis (m)',    'fontsize', font_size_label)
            ylabel('Y-axis (m)',    'fontsize', font_size_label)
            title( '3-DOF Robot',   'fontsize', font_size_title)

            n = 1;
            for(time=start_t:delta_t:finish_t)
                q1 = simul_q1(n);
                q2 = simul_q2(n);
                q3 = simul_q3(n);

                x1 = L1*cos(q1);               % [m]  : Joint 1 X-axis position
                y1 = L1*sin(q1);               % [m]  : Joint 1 Y-axis position
                x2 = L2*cos(q1+q2);       % [m]  : Joint 2 X-axis position
                y2 = L2*sin(q1+q2);       % [m]  : Joint 2 Y-axis position
                x3 = L3*cos(q1+q2+q3);       % [m]  : Joint 2 X-axis position
                y3 = L3*sin(q1+q2+q3);       % [m]  : Joint 2 Y-axis position

                Px1 = [0  x1];
                Py1 = [0  y1];
                Px2 = [x1 x1+x2];
                Py2 = [y1 y1+y2];
                Px3 = [x1+x2 x1+x2+x3];
                Py3 = [y1+y2 y1+y2+y3];

                set(p1, 'XData', Px1, 'YData', Py1)
                set(p2, 'XData', Px2, 'YData', Py2)
                set(p3, 'XData', Px3, 'YData', Py3)
                drawnow
                n = n+1;
            end
    end
    
    if(flag_Draw_Graph == 1)
        %Draw Postion
        FG2 = figure('Position',[900 700 600 300],'Color',[1 1 1]);
            % Target position visualize
            plot(simul_time,simul_q_d1 * 180/pi,':b','linewidth',linewidth_target); hold on;
            plot(simul_time,simul_q_d2 * 180/pi,':r','linewidth',linewidth_target); hold on;
            plot(simul_time,simul_q_d3 * 180/pi,':g','linewidth',linewidth_target); hold on;
            % Current position visualize
            plot(simul_time,simul_q1 * 180/pi,'b','linewidth',linewidth_current); hold on;
            plot(simul_time,simul_q2 * 180/pi,'r','linewidth',linewidth_current); hold on;
            plot(simul_time,simul_q3 * 180/pi,'g','linewidth',linewidth_current); hold on;
            
            axis([start_t finish_t -120 120]);
            xticks([start_t:0.5:finish_t])    % 눈금 그리는 것 (시작,눈금 간격, 끝)
            yticks([-120:30:120])              % 눈금 그리는 것 (시작,눈금 간격, 끝)
            grid on

            legend({'tar_{Link1}', 'tar_{Link2}','tar_{Link3}','cur_{Link1}','cur_{Link2}','cur_{Link3}'},'location','best','orientation','horizontal','fontsize',10)

            xlabel('time (s)',              'fontsize',font_size_label)
            ylabel('Angle (deg)',           'fontsize', font_size_label)
            title( 'Joint Space PID CTM Controller', 'fontsize', font_size_title)

       %Draw Angular Velocity

        FG3 = figure('Position',[900 300 600 300],'Color',[1 1 1]);
            % Target Angular Velocity visualize
            plot(simul_time,simul_dq_d1 * 180/pi,':b','linewidth',linewidth_target); hold on;
            plot(simul_time,simul_dq_d2 * 180/pi,':r','linewidth',linewidth_target); hold on;
            plot(simul_time,simul_dq_d3 * 180/pi,':g','linewidth',linewidth_target); hold on;
            % Current Angular Velocity visualize
            plot(simul_time,simul_dq1 * 180/pi,'b','linewidth',linewidth_current); hold on;
            plot(simul_time,simul_dq2 * 180/pi,'r','linewidth',linewidth_current); hold on;
            plot(simul_time,simul_dq3 * 180/pi,'g','linewidth',linewidth_current); hold on;

            axis([start_t finish_t -120 120]);
            xticks([start_t:0.5:finish_t])
            yticks([-120:20:120])
            grid on

            legend({'tar_{Link1}', 'tar_{Link2}','tar_{Link3}','cur_{Link1}','cur_{Link2}','cur_{Link3}'},'location','best','orientation','horizontal','fontsize',10)

            xlabel('time (s)',              'fontsize',font_size_label)
            ylabel('Angular Velocity (deg/s))',           'fontsize', font_size_label)
            title( 'Joint Space PID CTM Controller', 'fontsize', font_size_title)
    end
end
            
        
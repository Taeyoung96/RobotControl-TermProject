%Initialization
clc,clear all 
close all

%% Set Simulation Parameters    
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

        init_q1 = -pi/2; init_q2 = pi/4; init_q3 = pi/4;       %[rad]      : Init Joint Angle
        init_dq1 = 0.0;  init_dq2 = 0.0; init_dq3 = 0.0;       %[rad/s]    : Init Angular Velocity

        q           = [init_q1; init_q2; init_q3];      %[rad]      : Current Joint Angle
        dq          = [init_dq1; init_dq2; init_dq3];   %[rad/s]    : Current Angular Velocity

        pre_J = zeros(2,3);                              %          : pre jacobian initialization

        init_X  = get_Kinematics(q(1),q(2),q(3));     %[m]    : Init End-Effector Position
        
        X       = init_X;                        %[m]   : Current Position
        dX      = [0;0];                        %[m]    : Current Velocity
        
        X_d     = init_X;                        %[m]    : Current End-Effector Position
        dX_d    = [0;0];                        %[m/s]    : Current End-Effector Velocity
        
        ddX_d   = [0;0];                        %[m/s^2]    : Current End-Effector Acc

        pre_X = 0;
        dX_d_old = [0;0];    
        dX_d_int = [0;0];

        tq1 = 0.0;              % [Nm]  : Control Torque1
        tq2 = 0.0;              % [Nm]  : Control Torque2
        tq3 = 0.0;              % [Nm]  : Control Torque3
        tq = [tq1; tq2; tq3];              % [Nm]  : Control Torque
             
    %Controller Gain
        zeta        =   1;    %critical damped system
        Wn          =   20;         %[rad/s]        : natural frequency
        Kp          =   Wn^2;       %[Nm/rad]       : propotional Gain
        Kv          =   2*zeta*Wn;  %[Nm*s/rad]:    : Derivative Gain
        Ki          =   500;       %[Nm*s/rad]:    : Integrate Gain
            
%% Simulation
if(flag_Simul == 1)
    %Simulation
        n = 1;
        sin_t = 0;
        for(time = start_t:delta_t:finish_t)
        %Set Target Position  ??
            if(time < 1.0)
                X_d         = init_X;
                dX_d        = [0;0]; 
                ddX_d       = [0;0];
            elseif(time < 2.0)
                X_d(1)      = init_X(1);
                if(X_d(2) < init_X(2) + 0.1)
                    X_d(2) = X_d(2) + (0.1/0.5)*delta_t;
                else
                    X_d(2) = init_X(2) + 0.1;
                end
                dX_d = (X_d - [simul_X_d_x(n-1); simul_X_d_y(n-1)])./delta_t;
                ddX_d = (dX_d - [simul_dX_d_x(n-1); simul_dX_d_y(n-1)])./delta_t;
            else
                X_d     =[ 0.1 * sin((2*pi*sin_t)) + init_X(1);
                           0.1 * cos((2*pi*sin_t)) + init_X(2);];
                sin_t = sin_t + delta_t;
                dX_d    = (X_d - [simul_X_d_x(n-1); simul_X_d_y(n-1)])./delta_t;
                ddX_d   = (dX_d - [simul_dX_d_x(n-1); simul_dX_d_y(n-1)])./delta_t;
            end
                
                %Get Dynamics
                J = get_Jacobian(q(1),q(2),q(3));        %jacobian 
                %dJ = get_dJacobian(q(1),q(2));      %jacobian dot
                dJ = (J - pre_J)/delta_t;
                pre_J = J; 

                X = get_Kinematics(q(1),q(2),q(3));      % position
                dX = J*dq;                          % end effector velocity

                %dX = (X-pre_X)/delta_t;
                %pre_X = X;

                D = get_Inertia(q(2),q(3));                  %inertia
                H = get_Coriollis(q(2),q(3), dq(1),dq(2),dq(3));         %coriollis
                C = get_Gravity(q(1),q(2),q(3));             % gravity
                
                % For error function
                dX_d_int = dX_d_int + (dX_d - dX_d_old)*delta_t;
                dX_d_old = dX_d;

                %Control
                %u = ddX_d + Kp*(X_d - X);    % P control
                %u = ddX_d+Kv*(dX_d - dX) + Kp*(X_d - X);    % PD control
                u = ddX_d+Kv*(dX_d - dX) + Kp*(X_d - X) + Ki*dX_d_int;    % PID control
                ddq_ref = pinv(J)*(u - dJ*dq);           % 2 dot reference 

                tq_ctrl = D*ddq_ref + H + C*0.8;        % Torque control

             %Robot Model
                % Inverse Dynamics
                tq = tq_ctrl;
                tq1 = tq(1);    tq2 = tq(2);     tq3 = tq(3);
              
                [t,y] = ode45('three_links', [0 delta_t], [q(1); dq(1); q(2); dq(2); q(3); dq(3)]);
                index = length(y);
                q =  [y(index, 1); y(index,3); y(index,5)];
                dq = [y(index,2);  y(index,4); y(index,6)];

                %save Data
                simul_time(n)   = time;      %[sec]
                simul_q1(n)      = q(1);         %[rad]
                simul_q2(n)      = q(2);         %[rad]
                simul_q3(n)      = q(3);         %[rad]
                
                simul_dq1(n)     = dq(1);       % [rad/s]
                simul_dq2(n)      = dq(2);       % [rad/s]
                simul_dq3(n)      = dq(3);       % [rad/s]
                
                % Current
                simul_X_x(n)     = X(1);       % [m]
                simul_X_y(n)     = X(2);       % [m]

                simul_dX_x(n)     = dX(1);       % [m/s]
                simul_dX_y(n)     = dX(2);       % [m/s]
                
                % Desired
                simul_X_d_x(n)    = X_d(1);       % [m]
                simul_X_d_y(n)    = X_d(2);       % [m]

                simul_dX_d_x(n)   = dX_d(1);       % [m/s]
                simul_dX_d_y(n)   = dX_d(2);       % [m/s]

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

                axis([-1.0 1.0 -1.6 0.4]);
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

            plot(simul_time,simul_X_d_x,':r','linewidth',linewidth_target); hold on;
            plot(simul_time,simul_X_d_y,':b','linewidth',linewidth_target); hold on;

            plot(simul_time,simul_X_x,'r','linewidth',linewidth_current); hold on;
            plot(simul_time,simul_X_y,'b','linewidth',linewidth_current); hold on;

            axis([start_t finish_t -1.25 1]);
            xticks([start_t:1:finish_t])
            yticks([-1:0.25:1])
            grid on

            legend({'tar_x', 'tar_y','cur_x','cur_y'},'location','best','orientation','horizontal','fontsize',15)

            xlabel('time (s)',              'fontsize',font_size_label)
            ylabel('Position (m)',           'fontsize', font_size_label)
            title( 'Cartesin Space PID CTM Controller', 'fontsize', font_size_title)

        %Draw Velocity

        FG3 = figure('Position',[900 300 600 300],'Color',[1 1 1]);

            plot(simul_time,simul_dX_d_x,':r','linewidth',linewidth_target); hold on;
            plot(simul_time,simul_dX_d_y,':b','linewidth',linewidth_target); hold on;

            plot(simul_time,simul_dX_x,'r','linewidth',linewidth_current); hold on;
            plot(simul_time,simul_dX_y,'b','linewidth',linewidth_current); hold on;

            axis([start_t finish_t -1.25 1.25]);
            xticks([start_t:1:finish_t])
            yticks([-1.25:0.25:1.25])
            grid on

            legend({'tar_x', 'tar_y','cur_x','cur_y'},'location','best','orientation','horizontal','fontsize',15)

            xlabel('time (s)',              'fontsize',font_size_label)
            ylabel('Velocity (m/s)',           'fontsize', font_size_label)
            title( 'Cartesin Space PID CTM Controller', 'fontsize', font_size_title)

    end
end
        
        
        
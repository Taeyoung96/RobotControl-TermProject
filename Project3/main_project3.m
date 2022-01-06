%Initialization
clc,clear all
close all
%% Set Simulation Parameters      
%Draw flag
    flag_Simul       =1;
    flag_Draw        =1;
    flag_Draw_Robot = 1;
    flag_Draw_Graph = 1;

%Global Variable Initalization
    global I L G tq g m;

%simulation Parameters
    delta_t     = 0.005;    %[sec]      :Sampling Time
    start_t     = 0.000;    %[sec]      :Start Time
    finish_t    = 5.000;    %[sec]      : End Time

    g           =9.8148;    %[m/s^2]    : Gravitational Acceleration
    
    q_err_sum = 0;

    %Robot Parameters
    m           =1.0000;    %[kg]       :Link Inertia
    L           =1.0000;    %[m]       :Link Length
    I           =(m*L^2)/3; %[kgm^2]    :Link Inertia
    tq          =0.0000;    %[Nm]       :Control Torque

    init_q      =0;         %[rad]      : Init Joint Angle
    init_dq     =0;         %[rad/s]    : Init Angular Velocity

    q           =init_q;    %[rad]      : Current Joint Angle
    dq          = init_dq;   %[rad/s]    : Current Angular Velocity

% Target Position Parameters

    q_d         =   init_q;     %[rad]          : Target Joint Angle
    dq_d        =   init_dq;          %[rad/s]        : Target Angular Velocity
    ddq_d       =   0;          %[rad/s^2]      : Target Angular Acceleration

%Controller Gain
    zeta        = 1;    %critical damped system
    Wn          =   20;         %[rad/s]        : natural frequency
    Kp          =   Wn^2;       %[Nm/rad]       : propotional Gain
    Kv          =   2*zeta*Wn;  %[Nm*s/rad]:    : Derivative Gain
    Ki          =   2000;       %Nm*s/rad]:     : Integral Gain
            
%% Simulation
    if(flag_Simul == 1)
        %Simulation
            n = 1;
            for(time = start_t:delta_t:finish_t)
                %Set Target Trajectory
                    if(time < 1)
                        q_d         = init_q;
                        dq_d        = 0.0;
                        ddq_d       = 0.0;
                    else
                        if(q_d < 90*pi/180) % 현재 Desired 값이 90도 보다 작으면
                            q_d = q_d + (45*(pi/180)/2)*delta_t;      
                        else
                            q_d = 90*pi/180;
                        end
                        dq_d    = 30*pi/180;    %(q_d - simul_q_d(n-1))/delta_t;          %simul_q_d = pre_q_d
                        ddq_d   = (dq_d - simul_dq_d(n-1))/delta_t;        %simul_dq_d = pre_dq_d
                    end
                
                % Error Sum 
                q_err_sum = (q_d - q)*delta_t + q_err_sum;
                       
                % Get Dynamics   ---   1
                     G       = get_Gravity(q);         %get_Gravity(theta)
                % Controller
                  u        = ddq_d + Kv*(dq_d - dq) + Kp*(q_d - q) + Ki*q_err_sum; % PID controller
                  %u        = ddq_d + Kv*(dq_d - dq) + Kp*(q_d - q); %PD controller
                  %u        = ddq_d + Kp*(q_d - q);      %P controller
                                        
                  tq_ctrl = I*u + G*0.8; 
                %robot model
                    %Inverse Dynamics
                        tq      = tq_ctrl;
                        [t,y]   = ode45('one_link',[0 delta_t],[q; dq]);
                        index   = length(y);
                        q       =y(index,1);
                        dq      =y(index,2);
                    %save Data
                        simul_time(n)   = time;      %[sec]
                        simul_q(n)      = q;         %[rad]        joint angle
                        simul_dq(n)     = dq;       % [rad/s]
                        
                        simul_q_d(n)      = q_d;         %[rad]        
                        simul_dq_d(n)     = dq_d;       % [rad/s]
                        n               = n+ 1;
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
                init_x      = L*sin(init_q);
                init_y      = -L*cos(init_q);

                FG1 = figure('Position',[200 300 700 700], 'Color', [1 1 1]);
                    AX = axes('parent',FG1); hold on

                    p = plot([0 0], [init_x init_y], '-ob', 'Linewidth', linewidth_current);

                    axis([-1.5 1.5 -1.5 1.5]);
                    grid on
                xlabel('X-axis (m)',    'fontsize', font_size_label)
                ylabel('Y-axis (m)',    'fontsize', font_size_label)
                title( '1-DOF Robot',   'fontsize', font_size_title)
                
                n = 1;
                for(time=start_t:delta_t:finish_t)
                    q = simul_q(n);
                    x = L*sin(q); y = -L*cos(q);
                    Px = [0,x];     Py = [0,y];
                    set(p, 'XData', Px, 'YData', Py)
                    drawnow
                    n = n+1;
                end
        end

        if(flag_Draw_Graph == 1)
            %Draw Angle
                FG2 = figure('Position', [900 700 600 300], 'Color', [1 1 1]);
                plot(simul_time, simul_q_d * 180/pi,  ':k','linewidth', linewidth_target); hold on;
                plot(simul_time, simul_q * 180/pi,  ' r','linewidth', linewidth_current); hold on;

                legend('Desired','Current');
                axis([start_t finish_t 0 120]);
                xticks([start_t:1:finish_t])
                yticks([0:45:90])
                grid on
            xlabel('time (s)',              'fontsize',font_size_label)
            ylabel('Angle (deg)',           'fontsize', font_size_label)
            title( 'Joint Space PID CTM Controller', 'fontsize', font_size_title)

          %Draw Angular Velocity
            FG3 = figure('Position', [900 300 600 300], 'Color', [1 1 1]);
                plot(simul_time, simul_dq_d * 180/pi,     ':k', 'linewidth',linewidth_target); hold on;
                plot(simul_time, simul_dq * 180/pi,     '  r', 'linewidth',linewidth_current); hold on;

                legend('Desired','Current');
                axis([start_t finish_t -90 90]);
                xticks([start_t:1:finish_t])
                yticks([-60 30 60])
                grid on
            xlabel('time (s)',              'fontsize',font_size_label)
            ylabel('Angular Velocity (deg/s))',           'fontsize', font_size_label)
            title( 'Joint Space PID CTM Controller', 'fontsize', font_size_title)
        end
    end
        
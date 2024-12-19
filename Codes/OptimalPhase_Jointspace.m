% Solves the Optimization problem to guarantee joint velocities/accelerations
% remain bounded on specified limits
%
% To run this code you need first to run 'TrajParam_Jointspace.m'
% and download CasADi solver from https://web.casadi.org/
%
% credits for 'CasADi':
% 'CasADi: a software framework for nonlinear optimization and optimal control'
% DOI: 10.1007/s12532-018-0139-4
%
% Written by Giovanni Braglia, Davide Tebaldi and Luigi Biagiotti, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/


%--------------------------------------------------------------------------
function OptimalPhase_Jointspace
%--------------------------------------------------------------------------
    clear; clc; close all;

% PATH SHOULD BE CHANGED ACCORDING TO CasADi FOLDER
    addpath( '~/casadi-3.6.3-linux64-matlab2018b' );
    import casadi.*


    % Franka Joint Velocities/Accelerations Limits
    % https://frankaemika.github.io/docs/control_parameters.html
    dq_max  = [ 2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100 ]; % [rad/s]
    ddq_max = [ 15, 7.5, 10, 12.5, 15, 20, 20 ]; % [rad/s^2]

    addpath('Data');
    load('Data/Q_dim.mat');
    load('Data/Filtered_Trajectory_X.mat');
    xn = Q_dim;


    % Simulation Parameters
    param.nbData = length(xn);                             ... number of data samples to consider
    param.T      = tn(end);                                ... duration of demonstrated trajectory [s]
    param.Ts     = 0.001;                                  ... Sampling Time
    param.s0     = 0;                                      ... initial phase value
    param.sT     = sn(end);                                ... final phase value
    param.s      = [0; sn];                                ... phase values for numerical evaluations
    param.t      = linspace(0, param.T, param.nbData+1)';  ... time values for numerical evaluations
    param.DoF    = 7;                                      ... DoF of the robot


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opti = casadi.Opti(); % Optimization problem

    % ---- decision variables ---------
    X = opti.variable(3,param.nbData+1); % state trajectory
    s   = X(3,:);
    ds  = X(2,:);
    dds = X(1,:);
    U = opti.variable(1,param.nbData);   % control trajectory (JERK)
    T = opti.variable();                 % final time

    % ---- objective ------------------
    opti.minimize( T );

    % ---- dynamic constraints --------
    dt = T/param.nbData; % length of a control interval
    for k=1:param.nbData % loop over control intervals
    % chain of 3 intergrators
    Ad = [1 0 0; dt 1 0; 0 dt 1];
    Bd = [dt 0 0]';
    Cd = [0 0 1];
    opti.subject_to(X(:,k+1)==Ad*X(:,k) + Bd*U(:,k));
    end




    % ---- path constraints -----------
    opti.subject_to( -dq_max(1)  <= dq1_of_s(s).*ds  <= dq_max(1) );
    opti.subject_to( -ddq_max(1) <= ddq1_of_s(s).*ds.^2 + dq1_of_s(s).*dds  <= ddq_max(1) );
    opti.subject_to( -dq_max(2)  <= dq2_of_s(s).*ds <= dq_max(2) );
    opti.subject_to( -ddq_max(2) <= ddq2_of_s(s).*ds.^2 + dq2_of_s(s).*dds  <= ddq_max(2) );
    opti.subject_to( -dq_max(3)  <= dq3_of_s(s).*ds <= dq_max(3) );
    opti.subject_to( -ddq_max(3) <= ddq3_of_s(s).*ds.^2 + dq3_of_s(s).*dds  <= ddq_max(3) );
    opti.subject_to( -dq_max(4)  <= dq4_of_s(s).*ds <= dq_max(4) );
    opti.subject_to( -ddq_max(4) <= ddq4_of_s(s).*ds.^2 + dq4_of_s(s).*dds  <= ddq_max(4) );
    opti.subject_to( -dq_max(5)  <= dq5_of_s(s).*ds <= dq_max(5) );
    opti.subject_to( -ddq_max(5) <= ddq5_of_s(s).*ds.^2 + dq5_of_s(s).*dds  <= ddq_max(5) );
    opti.subject_to( -dq_max(6)  <= dq6_of_s(s).*ds <= dq_max(6) );
    opti.subject_to( -ddq_max(6) <= ddq6_of_s(s).*ds.^2 + dq6_of_s(s).*dds  <= ddq_max(6) );
    opti.subject_to( -dq_max(7)  <= dq7_of_s(s).*ds <= dq_max(7) );
    opti.subject_to( -ddq_max(7) <= ddq7_of_s(s).*ds.^2 + dq7_of_s(s).*dds <= ddq_max(7) );
    opti.subject_to( -100 <= U <= 100 ); % control is limited


    % ---- boundary conditions --------
    opti.subject_to( s(1) == param.s0 );
    opti.subject_to( ds(1) == 0 );
    opti.subject_to( dds(1) == 0 );
    opti.subject_to( s(param.nbData+1) == param.sT );
    opti.subject_to( ds(param.nbData+1) == 0 );
    opti.subject_to( dds(param.nbData+1) == 0 );

    % ---- misc. constraints  ----------
    opti.subject_to( T > 0 ); % Time must be positive
    opti.subject_to( T < param.T );
    opti.subject_to( ds >= 0 );
    ...opti.subject_to( param.s0 <= s <= param.sT );


    % ---- initial values for solver ---
    opti.set_initial( s, param.s0 );
    opti.set_initial( ds, 0 );
    opti.set_initial( dds, 0 );
    opti.set_initial( T, 1 );

    % ---- solve NLP              ------
    opti.solver('ipopt');
    sol = opti.solve();   % actual solve

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    T_opt   = sol.value(T);
    s_opt   = sol.value(s);
    ds_opt  = sol.value(ds);
    dds_opt = sol.value(dds);

    s_opt_ = s_opt;

    Q_opt = [                 ...
            q1_of_s(s_opt_)' ,...
            q2_of_s(s_opt_)' ,...
            q3_of_s(s_opt_)' ,...
            q4_of_s(s_opt_)' ,...
            q5_of_s(s_opt_)' ,...
            q6_of_s(s_opt_)' ,...
            q7_of_s(s_opt_)' ,...
            ];

    % Create arrays for final plots
    Q = [                   ...
            q1_of_s(param.s) ,...
            q2_of_s(param.s) ,...
            q3_of_s(param.s) ,...
            q4_of_s(param.s) ,...
            q5_of_s(param.s) ,...
            q6_of_s(param.s) ,...
            q7_of_s(param.s) ,...
        ];

    dQ_opt = [                 ...
            dq1_of_s(s_opt_)' ,...
            dq2_of_s(s_opt_)' ,...
            dq3_of_s(s_opt_)' ,...
            dq4_of_s(s_opt_)' ,...
            dq5_of_s(s_opt_)' ,...
            dq6_of_s(s_opt_)' ,...
            dq7_of_s(s_opt_)' ,...
            ];

    ddQ_opt = [                 ...
            ddq1_of_s(s_opt_)' ,...
            ddq2_of_s(s_opt_)' ,...
            ddq3_of_s(s_opt_)' ,...
            ddq4_of_s(s_opt_)' ,...
            ddq5_of_s(s_opt_)' ,...
            ddq6_of_s(s_opt_)' ,...
            ddq7_of_s(s_opt_)' ,...
            ];


    time = linspace( 0, T_opt, param.nbData+1 );
    tf   = linspace( T_opt, tn(end), param.nbData );
    sf   = sn(end) .* ones( [1,param.nbData] );

    figure(1); clf;
    %-----
    subplot(3,1,1)
    plot( time, s_opt, 'LineWidth', 2.5 );
    xlabel('time [s]')
    ylabel('s_{opt}')

    subplot(3,1,2)
    plot( time, ds_opt, 'LineWidth', 2.5 );
    xlabel('time [s]')
    ylabel('ds_{opt}')

    subplot(3,1,3)
    plot( time, dds_opt, 'LineWidth', 2.5 );
    xlabel('time [s]')
    ylabel('dds_{opt}')


    figure(2); clf;
    %-----
    plot( tn, sn, 'LineWidth', 2, 'DisplayName', 's demonstration' );
    hold on;
    plot( [time tf], [s_opt_ sf], 'LineWidth', 2, 'DisplayName', 's optimal' );
    xlabel('time [s]')
    ylabel('s')
    leg = legend();
    set(leg,'FontSize',10);

    figure(3); clf;
    %-----
    subplot(2,1,1)
    plot( param.t, Q(:,1), 'LineWidth', 2 );
    hold on
    grid on
    for i=2:param.DoF
        plot( param.t, Q(:,i), 'LineWidth', 2 );
    end
    xlabel('time [s]')
    ylabel('$\textbf{q}$ [rad]', 'Interpreter', 'latex')
    leg = legend('$q_1$', '$q_2$', '$q_3$', '$q_4$', '$q_5$', '$q_6$', '$q_7$');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);

    subplot(2,1,2)
    plot( time, Q_opt(:,1), 'LineWidth', 2 );
    hold on
    grid on
    for i=2:param.DoF
        plot( time, Q_opt(:,i), 'LineWidth', 2 );
    end
    xlabel('time [s]')
    ylabel('$\hat{\textbf{q}}$ [rad]', 'Interpreter', 'latex')
    leg = legend('$\hat{q}_1$', '$\hat{q}_2$', '$\hat{q}_3$', '$\hat{q}_4$', '$\hat{q}_5$', '$\hat{q}_6$', '$\hat{q}_7$' );
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);


    figure(4); clf;
    %-----
    subplot(2,1,1)
    plot( time, (dQ_opt(:,1).*ds_opt')/dq_max(1) , 'LineWidth', 2.5 );
    hold on
    grid on
    for i=2:param.DoF
        plot( time, (dQ_opt(:,i).*ds_opt')/dq_max(i), 'LineWidth', 2.5 );
    end
    xlabel('time [s]')
    ylabel('Normalized $\hat{\dot{\textbf{q}}}$ [$rad/s$]', 'Interpreter', 'latex')
    leg = legend('$\hat{\dot{q}}_1$', '$\hat{\dot{q}}_2$', '$\hat{\dot{q}}_3$', '$\hat{\dot{q}}_4$', '$\hat{\dot{q}}_5$', '$\hat{\dot{q}}_6$', '$\hat{\dot{q}}_7$' );
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);

    subplot(2,1,2)
    plot( time, (ddQ_opt(:,1).*ds_opt'.^2 + dQ_opt(:,1).*dds_opt')/ddq_max(1) , 'LineWidth', 2.5 );
    hold on
    grid on
    for i=2:param.DoF
        plot( time, (ddQ_opt(:,i).*ds_opt'.^2 + dQ_opt(:,i).*dds_opt')/ddq_max(i), 'LineWidth', 2.5 );
    end
    xlabel('time [s]')
    ylabel('Normalized $\hat{\ddot{\textbf{q}}}$ [$rad/s^2$]', 'Interpreter', 'latex')
    leg = legend('$\hat{\ddot{q}}_1$', '$\hat{\ddot{q}}_2$', '$\hat{\ddot{q}}_3$', '$\hat{\ddot{q}}_4$', '$\hat{\ddot{q}}_5$', '$\hat{\ddot{q}}_6$', '$\hat{\ddot{q}}_7$' );
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);

end

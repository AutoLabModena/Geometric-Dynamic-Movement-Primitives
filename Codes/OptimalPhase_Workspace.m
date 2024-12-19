% Solves the Optimization problem to guarantee constant feedrate
% on Cartesian velocity
%
% To run this code you need first to run 'TrajParam_Workspace.m'
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
function OptimalPhase_Workspace
%--------------------------------------------------------------------------
clear; clc;

% PATH SHOULD BE CHANGED ACCORDING TO CasADi FOLDER
    addpath( '~/casadi-3.6.3-linux64-matlab2018b' );
    import casadi.*

    addpath('Data');
    load('Data/Filtered_Trajectory_X.mat'); % from 'TrajParam_Workspace.m'

    % Franka WorkSpace Velocities/Accelerations Limits
    vmax = 0.5;  % [m/s]
    amax = 10.0; % [m/s^2]


    % Simulation Parameters
    param.nbData = length(xn);                                       ... number of data samples to consider
    param.T      = tn(end);                                          ... duration of demonstrated trajectory [s]
    param.Ts     = 0.001;                                            ... Sampling Time
    param.s0     = 0;                                                ... initial phase value
    param.sT     = sn(end);                                          ... final phase value
    param.s      = [0; sn];                                          ... phase values for numerical evaluations
    param.t      = linspace(0, param.T, param.nbData+1)';            ... time values for numerical evaluations
    param.DoF    = 3;                                                ... DoF of the robot


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opti = casadi.Opti(); % Optimization problem

    % ---- decision variables --------- +

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

    % delta = 0.05;
    % ---- path constraints -----------
    opti.subject_to( -vmax < ds  < vmax );
    opti.subject_to( -amax  < dds < amax );
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
    opti.subject_to( ds >= 0 );
    opti.subject_to( T < param.T );
    opti.subject_to( param.s0 <= s <= param.sT );


    % ---- initial values for solver ---
    opti.set_initial( s, param.s0 );
    opti.set_initial( ds, 0 );
    opti.set_initial( dds, 0 );
    opti.set_initial( T, 1 );

    % ---- solve NLP              ------
    opti.solver('ipopt'); % set numerical backend
    sol = opti.solve();   % actual solve

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    T_opt   = sol.value(T);
    s_opt   = sol.value(s);
    ds_opt  = sol.value(ds);
    dds_opt = sol.value(dds);


    s_opt_ = s_opt;

    P_opt  = [                 ...
                x1_of_s( s_opt )', ...
                x2_of_s( s_opt )', ...
                x3_of_s( s_opt )'  ...
            ];

    P  = [             ...
            x1_of_s( param.s ), ...
            x2_of_s( param.s ), ...
            x3_of_s( param.s )  ...
        ];


    V_opt  = [                  ...
                dx1_of_s( s_opt )', ...
                dx2_of_s( s_opt )', ...
                dx3_of_s( s_opt )'  ...
            ];

    V  = [              ...
            dx1_of_s( param.s ), ...
            dx2_of_s( param.s ), ...
            dx3_of_s( param.s )  ...
        ];



    time = linspace( 0, T_opt, param.nbData+1 );
    tf   = linspace( T_opt, tn(end), param.nbData );
    sf   = sn(end) .* ones( [1,param.nbData] );

    figure
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

    figure
    %-----
    plot( tn, sn, 'LineWidth', 2.5, 'DisplayName', 's demonstration' );
    hold on;
    plot( [time tf], [s_opt_ sf], 'LineWidth', 2.5, 'DisplayName', 'constant feedrate' );
    xlabel('time [s]')
    ylabel('s')
    leg = legend();
    set(leg,'FontSize',10);

    figure
    %-----
    subplot(2,1,1)
    plot( param.t, P(:,1), 'LineWidth', 2.5 );
    hold on
    grid on
    for i=2:param.DoF
        plot( param.t, P(:,i), 'LineWidth', 2.5 );
    end
    xlabel('time [s]')
    ylabel('$\textbf{p}$ [m]', 'Interpreter', 'latex')
    leg = legend('$x$', '$y$', '$z$');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);

    subplot(2,1,2)
    plot( time, P_opt(:,1), 'LineWidth', 2.5 );
    hold on
    grid on
    for i=2:param.DoF
        plot( time, P_opt(:,i), 'LineWidth', 2.5 );
    end
    xlabel('time [s]')
    ylabel('$\hat{\textbf{p}}$ [m]', 'Interpreter', 'latex')
    leg = legend('$\hat{x}$', '$\hat{y}$', '$\hat{z}$');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);


    normV = zeros( [ length(V_opt), 1 ] );
    vM    = ones( [ length(V_opt), 1 ] );
    for i=1:length(V)
        normV(i) = (norm( V_opt(i,:) ) * ds_opt(i))/vmax;
    %     vM(i)    = vmax;
    end
    figure
    %-----

    plot( time, abs(ds_opt./vmax), 'b', 'LineWidth', 2.5 );
    hold on
    plot( time, vM,'r--', 'LineWidth', 2.5 );
    xlabel('time [s]')
    % ylabel('$|| \hat{\dot{\textbf{p}}} ||$', 'Interpreter', 'latex')
    leg = legend( 'Normalized $|| \hat{\dot{\textbf{p}}} ||$', 'Normalized $v_{MAX}$');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);

    figure
    %-----
    ddsM  = ones( [ length(dds_opt), 1 ] );
    plot( time, dds_opt./amax, 'b', 'LineWidth', 2.5 );
    hold on
    plot( time, ddsM,'r--', 'LineWidth', 2.5 );
    xlabel('time [s]')
    % ylabel('$|| \hat{\dot{\textbf{p}}} ||$', 'Interpreter', 'latex')
    leg = legend( 'Normalized $\hat{\ddot{s}}$', 'Normalized $\ddot{s}_{MAX}$');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',12);

end

%
% Creates the joint trajectories 'Q_dim.mat' by inverse kinematics
% of the spatial-filtered workspace trajectories
%
% Written by Giovanni Braglia, Davide Tebaldi and Luigi Biagiotti, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/
%
function Workspace2Jointspace
    
    clear; clc;

    addpath('Data');
    load('Data/X_dim.mat');
    
    l     = length(X_dim);
    Ts    = 0.001; % Sampling time of recorded trajectory
    
    ptr      = X_dim;
    time_ptr = linspace(0,l*Ts,l)'; 
    delta    = 0.006;
    
    % Spatial Sampling of the recorded trajectory
    [tn,~,xn]= SpatialSampling( time_ptr, ptr, delta );
    
    
    xtr = xn(:,1);
    ytr = xn(:,2);
    ztr = xn(:,3);
    
    Q_dim = zeros(length(xtr),7);
    
    panda = importrobot('Data/panda_model.urdf'); 
    panda.Gravity = [0,0,-9.81];
    panda.Bodies{1,1}.Joint.HomePosition = 0;       % panda_joint1
    panda.Bodies{1,2}.Joint.HomePosition = -pi/4;   % panda_joint2
    panda.Bodies{1,3}.Joint.HomePosition = 0;       % panda_joint3
    panda.Bodies{1,4}.Joint.HomePosition = -3*pi/4; % panda_joint4
    panda.Bodies{1,5}.Joint.HomePosition = 0;       % panda_joint5
    panda.Bodies{1,6}.Joint.HomePosition = pi/2;    % panda_joint6
    panda.Bodies{1,7}.Joint.HomePosition = pi/4;    % panda_joint7
    
    params.MaxIterations = 250;
    params.StepTolerance = 10^-7;
    
    % Impose the homeConfiguration:
    panda.DataFormat = 'row';
    q_home = homeConfiguration(panda);
    
    % Create the System object for solving generalized inverse kinematics:
    gik = generalizedInverseKinematics;
    gik.SolverParameters = params;
    gik.RigidBodyTree = panda;
    
    gik.ConstraintInputs = {'pose'};
    poseTgt = constraintPoseTarget('panda_hand');
    
    for ii=1:length(xtr)
    
     pos =[xtr(ii)    ytr(ii)    ztr(ii)];
     poseTgt.TargetTransform = [ 1  0  0  pos(1);
                                 0 -1  0  pos(2);
                                 0  0 -1  pos(3);
                                 0  0  0   1];
    
     % Find solution:
     [qRobotics,~] = gik(q_home,poseTgt);
    
     % Search solution that doesn't reach the joint limits:
     Q_dim(ii,:) = qRobotics(1:7);
     q_home      =  q_home*0;
     q_home(1:7) = qRobotics(1:7);
     q_home(8:9) = [0 0];
    
     clc;
    
    end
    
    qn = Q_dim;

    figure;
    subplot(2,1,1)
    plot( tn, xn, 'LineWidth', 2 );
    xlabel('Time [s]')
    ylabel('m')
    title('Workspace trajectories, spatial filtered');

    subplot(2,1,2)
    plot( tn, qn, 'LineWidth', 2 );
    xlabel('Time [s]')
    ylabel('rad')
    title('Jointspace trajectories, spatial filtered');

    save('Data/Q_dim.mat', 'Q_dim');
end

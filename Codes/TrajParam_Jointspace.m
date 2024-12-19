% Calculates parameterized symbolic function of the recorded trajectory 'Q_dim.mat'.
%
% The recorded trajectory is a [DoF+1, N] array, the first row
% is the recorded time (not used), N is the number of datapoints.
%
% To run this code first run 'TrajParam_Workspace.m' then
% 'Workspace2Jointspace.m'
%
% Written by Giovanni Braglia and Davide Tebaldi, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/

clear all
close all
clc


DoF = 7; % Degree of Freedom, 7 for Franka robot

addpath('Data');
load('Data/Q_dim.mat');
load('Data/Filtered_Trajectory_X.mat');

qn = Q_dim;
Ts = 0.001; % Sampling time



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NrFig=1;
Fig_Name='Figure_1';
hh=figure(NrFig); clf;
set(hh,'NumberTitle','off')
set(hh,'Name',[num2str(NrFig) '.' Fig_Name])
% set(hh,'PaperPositionMode','manual','PaperPosition',[5 0 15 17.75]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:DoF
    subplot(4,2,ii)
    grid on; zoom on; hold on;
    plot(tn,qn(:,ii),'LineWidth',2)
    xlabel('Time [s]')
    ylabel('rad')
    title(['q' num2str(ii)])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NrFig=2;
Fig_Name='Figure_2';
hh=figure(NrFig); clf;
set(hh,'NumberTitle','off')
set(hh,'Name',[num2str(NrFig) '.' Fig_Name])
% set(hh,'PaperPositionMode','manual','PaperPosition',[5 0 15 17.75]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_range=sn;
s_fin=s_range(end);
time_qtr=tn;
ds=diff(s_range)/(time_qtr(2)-time_qtr(1));
ds2=ds.^2;
dds=diff(ds)/(time_qtr(2)-time_qtr(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameterization with Radial Basis Functions (RBF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 22;  % Number of RBF
L = length( qn );
c = linspace ( 0 - 1/(N-3), s_fin + 1/(N-3), N )' ; % centers of RBF

h = zeros( size(c) );
k = 1.3; ... [1.2, 1.5]
for i=1:N-1
    h(i) = 0.5* ( k * ( c(i+1) - c(i) ) )^2;
end
h(end) = h(end-1);


RBF   = cell ( [ N, 1] );    ... array of RBF functions
dRBF  = cell ( [ N, 1] );    ... array of first derivatives    
ddRBF = cell ( [ N, 1] );    ... array of second derivatives

for i = 1 : N 
    RBF{i}   = @(s) exp( - ( s - c(i) ).^2 / ( 2 * h(i) ) ); 
    dRBF{i}  = @(s) ( - ( s - c(i) ) / h(i) ) .* exp( - ( s - c(i) ).^2 / ( 2 * h(i) ) );
    ddRBF{i} = @(s) ( -1/h(i) + ( -( s - c(i) ) / h(i) ) .^2 ) .* exp( - ( s - c(i) ).^2 / ( 2 * h(i) ) );
end    

sum_rbf   = @(s) RBF{1}(s);
dsum_rbf  = @(s) dRBF{1}(s);
ddsum_rbf = @(s) ddRBF{1}(s);
for i = 2 : N 
    sum_rbf   = @(s) sum_rbf(s)   + RBF{i}(s);
    dsum_rbf  = @(s) dsum_rbf(s)  + dRBF{i}(s);
    ddsum_rbf = @(s) ddsum_rbf(s) + ddRBF{i}(s);
end

Phi   = cell ( [ N, 1] ); 
dPhi  = cell ( [ N, 1] ); 
ddPhi = cell ( [ N, 1] ); 

for i = 1:N
    Phi{i}   = @(s) RBF{i}(s) ./ sum_rbf(s);
    dPhi{i}  = first_der( RBF{i}, dRBF{i}, sum_rbf, dsum_rbf );
    ddPhi{i} = second_der( RBF{i}, dRBF{i}, ddRBF{i}, sum_rbf, dsum_rbf, ddsum_rbf );
end

qtr_par   = ones( [L,DoF] );
dqtr_par  = ones( [L,DoF] );
ddqtr_par = ones( [L,DoF] );

for ii=1:DoF

    qd=qn(:,ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PSI=ones([L,N]);
    dPSI=ones([L,N]);
    ddPSI=ones([L,N]);

    Ayg=ones([2,N]);   % constraint matrix for initial/final positions
    byg=[qd(1),qd(end)]';

    for jj=1:N
        PSI(:,jj)=Phi{jj}(s_range);
        dPSI(:,jj)=dPhi{jj}(s_range);
        ddPSI(:,jj)=ddPhi{jj}(s_range);
        %
        Ayg(1,jj)=Phi{jj}(0);
        Ayg(2,jj)=Phi{jj}(s_fin);
    end


    dqd  = diff(qd)/(Ts);
    ddqd = diff(dqd)/(Ts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weights Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W = LagrangeApprox( qd, PSI, Ayg, byg, N );
    clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create symbolic functions for Simulink ( NORMALIZED VERSION )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%   | ---------------------------------------- |
%   | The followings creates symbolic          |
%   | 'file.m' to be used in Simulink or       |
%   | to test  'OptimalPhase_Jointspace.m'.    |
%   | ---------------------------------------- |
%
    q   = @(s) W(1).*Phi{1}(s);
    dq  = @(s) W(1).*dPhi{1}(s);
    ddq = @(s) W(1).*ddPhi{1}(s);

    for jj = 2:N
        q   = @(s) q(s)   + W(jj).*Phi{jj}(s);
        dq  = @(s) dq(s)  + W(jj).*dPhi{jj}(s);
        ddq = @(s) ddq(s) + W(jj).*ddPhi{jj}(s);
    end

    ysym   = q(sym('s'));
    matlabFunction(ysym,'File',['Data/q' num2str(ii) '_of_s']);

    dysym  = dq(sym('s'));
    matlabFunction(dysym,'File',['Data/dq' num2str(ii) '_of_s']);

    ddysym = ddq(sym('s'));
    matlabFunction(ddysym,'File',['Data/ddq' num2str(ii) '_of_s']);

    clear q dq ddq ysym dysym ddysym


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    qtr_par(:,ii)=PSI*W;
    dqtr_par(:,ii)=dPSI*W;
    ddqtr_par(:,ii)=ddPSI*W;

    figure(2)
    %
    subplot(3,7,ii);
    plot(s_range,qtr_par(:,ii),'b' ,'LineWidth',2);
    hold on; grid on; zoom on;
    plot(s_range,qd,'r--','LineWidth',2);
    ylabel('rad')
    title('actual (r) and par (b)')
    %
    subplot(3,7,ii+7);
    hold on; grid on; zoom on;
    plot(s_range,dqtr_par(:,ii),'b','LineWidth',2);
    ylabel('rad/s')
    title('par vel')
    %
    subplot(3,7,ii+14);
    hold on; grid on; zoom on;
    plot(s_range,ddqtr_par(:,ii),'b','LineWidth',2);
    xlabel('Phase s')
    ylabel('rad/s^2')
    title('par acc')


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save variables for Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Data/Filtered_Trajectory_Q.mat','tn','sn','qn')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------------------------
function [tn,sn,xn]= SpatialSampling(t,x,delta)
%---------------------------------------------------------------------
% Returns equally spaced ('delta') points on trajectory 'x'

    [nrow,ncol]=size(x);
    if ncol>nrow
        x=x';
    end
    %
    for jj=1:size(x,2)
        eval(['xn' num2str(jj) '=[x(1,' num2str(jj) ')];']);
    end
    %
    sn=[t(1)];
    tn=[t(1)];
    %
    xcurr=x(1,:);
    for ii=1:size(x,1)-1
        %
        if norm(x(ii+1,:)-xcurr)>delta
            nr_steps=floor(norm(x(ii+1,:)-xcurr)/delta);
            for pp=1:nr_steps
                sn=[sn; sn(end)+delta];
                for jj=1:size(x,2)
                    xincr=xcurr(jj)+delta*(x(ii+1,jj)-xcurr(jj))/norm(x(ii+1,:)-xcurr);
                    eval(['xn' num2str(jj) '=[xn' num2str(jj) '; ' num2str(xincr) '];']);
                end
                xcurr=x(1,:)*0;
                for jj=1:size(x,2)
                    eval(['xcurr(' num2str(jj) ')=xn' num2str(jj) '(end);']);
                end
                nDX=norm(x(ii+1,:)-x(ii,:));
                DT=t(ii+1)-t(ii);
                K=nDX/DT;
                tn=[tn; (norm(xcurr-x(ii,:))+K*t(ii))/K];
            end
        end
    end
    xn=zeros(size(xn1,1),size(x,2));
    for jj=1:size(x,2)
        eval(['xn(:,' num2str(jj) ')=xn' num2str(jj) ';']);
    end

end

%---------------------------------------------------------------------
function dphi = first_der( rbf, drbf, sum_rbf, dsum_rbf )
%---------------------------------------------------------------------
% Calculates analitically the first derivative of normalized RBF

    num = @(s) drbf(s) .* sum_rbf(s) - rbf(s) .* dsum_rbf(s);
    den = @(s) sum_rbf(s).^2;

    dphi = @(s) num(s) ./ den(s);

end


%---------------------------------------------------------------------
function ddphi = second_der( rbf, drbf, ddrbf, sum_rbf, dsum_rbf, ddsum_rbf )
%---------------------------------------------------------------------
% Calculates analitically the second derivative of normalized RBF


    a = @(s) ddrbf(s) .* sum_rbf(s) - rbf(s) .* ddsum_rbf(s);
    b = @(s) drbf(s) .* sum_rbf(s) - rbf(s) .* dsum_rbf(s);


    num = @(s) a(s) .* ( sum_rbf(s).^2 ) -  b(s) .* ( 2 * sum_rbf(s) ) .* dsum_rbf(s);
    den = @(s) ( sum_rbf(s).^4 );

    ddphi = @(s) num(s) ./ den(s);

end


%---------------------------------------------------------------------
function P = LagrangeApprox( yd, PSI, Ayg, byg, N )
%---------------------------------------------------------------------
% Computes analitically the weighting coefficients for the
% approximation of 'yd' with RBF basis functions

    L  = length( yd );
    W  = eye( L ) .* 1e0;
    Bi = Ayg; 
    Ri = byg; 
    Ba = PSI;
    Ra = yd;

    Ba2 = Ba' * W * Ba + eye(N).*realmin;
    BRa = Ba' * W * Ra;

    [ k, ~ ] = size( Bi );
    Q = ( Bi / Ba2 ) * Bi'  + eye(k).*realmin ;

    l = Q \ ( ( Bi / Ba2 ) * BRa - Ri );
    P = ( Ba2 \ BRa ) - ( Ba2 \ Bi' ) * l ;

end

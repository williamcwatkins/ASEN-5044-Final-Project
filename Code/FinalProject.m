%% Admin
%{

Names: Corey LePine and William Watkins 
Professor: McMahon
Class: ASEN 5044 Stat Est for Dyn Sys
Date: December 14, 2021
Final Project: Statistical Orbit Determination

%} 
%% Housekeeping
clc; clear all; close all;

%% Constants
load('orbitdeterm_finalproj_KFdata.mat')
u = 398600; % Earth's standard gravitational paremters [km^3/s^2]
r0 = 6678; % Nominal orbit radius [km]
Re = 6378; % Uniform radius of Earth [km]
we = 2*pi/86400; % Consant rotation rate of Earth [rad/s]
vel = sqrt(u/r0); % orbital velocity
circ = 2*pi*r0; % circumference of the orbit
T = circ / vel; % Orbital period

n = 4; % number of states
m = 2; % number of inputs, also happens to be number of distrubances
p = 3; % number of measurements
j = 6; % number of measurement stations

dt = 10; % step size [s]
tspan = 0:dt:14000;
initCon = [r0, 0, 0, r0*sqrt(u/r0^3)]; % initial state point - LTV sys, so
% have to linearize about a nominal trajectory, see below

DEBUG = 0; % Debug boolean

ColorSet = varycolor(12); % 12 Unique Colors for the Tracking Stations

%% Part 1

% 1a) Find the CT Jacobian Matracies
syms x1 x2 x3 x4 mu u1 u2 w1 w2 z1 z2 z3 z4 t % zi is the states of the ground stations

f = [x2;
    -mu*x1/sqrt(x1^2 + x3^2)^3 + u1 + w1;
    x4;
    -mu*x3/sqrt(x1^2 + x3^2)^3 + u2 + w2];
state = [x1, x2, x3, x4];
inputs = [u1, u2];
disturb = [w1, w2];

A = jacobian(f, state);
B = jacobian(f, inputs);
Gam = jacobian(f, disturb);

h = [sqrt((x1 - z1)^2 + (x3 - z3)^2);
    ((x1 - z1)*(x2 - z2) + (x3 - z3)*(x4 - z4))/(sqrt((x1 - z1)^2 + (x3 - z3)^2));
    atan((x3 - z3)/(x1 - z1))];

C = jacobian(h, state); 
D = jacobian(h, inputs); 

% 1b) Linearize about Nominal operating points
for ii = 1:length(tspan)
    nomCon(:,ii) = [r0 * cos(sqrt(u/r0^3)*tspan(ii)); -r0 * sin(sqrt(u/r0^3)*tspan(ii))*sqrt(u/r0^3);
            r0*sin(sqrt(u/r0^3)*tspan(ii)); r0*cos(sqrt(u/r0^3)*tspan(ii))*sqrt(u/r0^3)];
            % Have to linearize about nominal trajectory!
end
mu = u;
    
%% 1c) Nonlinear Dynamics
% State NL 
Rel_Tol = 1e-13;
Abs_Tol = Rel_Tol;
options = odeset('Stats', 'off', 'RelTol', Rel_Tol, 'AbsTol', Abs_Tol);

perts = [0, 0.075, 0, -0.021];
Initial_States = perts + initCon;

[Time_out, State_out] = ode45(@(Time, State) StatODNL_ODE(Time, State), tspan, Initial_States, options);

State_X = State_out(:, 1);
State_Xdot = State_out(:, 2);
State_Y = State_out(:, 3);
State_Ydot = State_out(:, 4);

figure()
subplot(4, 1, 1)
plot(Time_out, State_out(:, 1))
xlabel('Time [s]')
ylabel('X [km]')
ylim([-1e4, 1e4])

subplot(4, 1, 2)
plot(Time_out, State_out(:, 2))
xlabel('Time [s]')
ylabel('Xdot [km/s]')

subplot(4, 1, 3)
plot(Time_out, State_out(:, 3))
xlabel('Time [s]')
ylabel('Y [km]')
ylim([-1e4, 1e4])

subplot(4, 1, 4)
plot(Time_out, State_out(:, 4))
xlabel('Time [s]')
ylabel('Ydot [km/s]')

sgtitle('States vs Time, Full Nonlinear Dynamics Simulation')

% Tracking Stations Positions
TS_IDS = 1:1:12; % tracking stations ids
theta_TS0 = (TS_IDS - 1)*pi/6; % tracking stations intial positions

for ii = 1:12 % tracking station X and Y position 
    TS_X(:, ii) = Re*cos(we*Time_out + theta_TS0(ii));
    TS_Xdot(:, ii) = -Re*we*sin(we*Time_out + theta_TS0(ii));
    TS_Y(:, ii) = Re*sin(we*Time_out + theta_TS0(ii));
    TS_Ydot(:, ii) = Re*we*cos(we*Time_out + theta_TS0(ii));
    theta_TS(:, ii) = atan2(TS_Y(:, ii), TS_X(:, ii));
end

thetaCompare = theta_TS;

% Measurement NL
for ii = 1:12
    rho(:, ii) = sqrt((State_X - TS_X(:, ii)).^2 + (State_Y - TS_Y(:, ii)).^2);
    rho_dot(:, ii) = ((State_X - TS_X(:, ii)).*(State_Xdot - TS_Xdot(:, ii)) + (State_Y - TS_Y(:, ii)).*(State_Ydot - TS_Ydot(:, ii)))./rho(:, ii);
    phi(:, ii) = atan2((State_Y - TS_Y(:, ii)), (State_X - TS_X(:, ii)));
    visibleStation(:,ii) = ones(length(tspan),1) * ii;
end

phiCompare = phi;

% Wrap the upper and lower bounds between -pi and pi
% When the upper bound is above pi, need to wrap both bound down to -pi
% Vice versa when lower bound is 
thetaBound1Pos = theta_TS;
thetaBound1PosInd = find(thetaBound1Pos+pi/2 > pi);
thetaBound1Pos(thetaBound1PosInd) = thetaBound1Pos(thetaBound1PosInd) - 2*pi;
thetaBound1Neg = theta_TS;
thetaBound1Neg(thetaBound1PosInd) = thetaBound1Neg(thetaBound1PosInd) - 2*pi;

thetaBound2Neg = theta_TS;
thetaBound2NegInd = find(thetaBound2Neg-pi/2 < pi);
thetaBound2Neg(thetaBound2NegInd) = thetaBound2Neg(thetaBound2NegInd) + 2*pi;
thetaBound2Pos = theta_TS;
thetaBound2Pos(thetaBound2NegInd) = thetaBound2Pos(thetaBound2NegInd) + 2*pi;

% Visible Tracking Stations
figure
hold on
tl = tiledlayout(4,1);
title(tl, "Full Nonlinear Model Data Simulation");
xlabel(tl, "Time (secs)");
nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    if(DEBUG == 1)
        plot(Time_out, pi/2 + thetaBound1Pos(:,ii));
        hold on
        plot(Time_out, -pi/2 + thetaBound1Neg(:,ii));
        plot(Time_out, phi(:,ii));
        plot(Time_out, theta_TS(:,ii));
        scatter(Time_out(vis_index), phi(vis_index, ii)) 
        yline(pi);
        yline(-pi);
    end
    scatter(Time_out(vis_index), rho(vis_index,ii));
    ylabel('rho^i (km)');
end

nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), rho_dot(vis_index,ii));
    ylabel('rhodot^i (km/s)');
end
nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), phi(vis_index,ii));
    ylabel('\phi^i (rads)');
end
nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), visibleStation(vis_index,ii));
    ylabel('Visible Station ID');
end


%% 1c) Linearized DT
pertX(:, 1) = perts';
for ii = 2:1401
    x1 = nomCon(1,ii-1);
    x2 = nomCon(2,ii-1);
    x3 = nomCon(3,ii-1);
    x4 = nomCon(4,ii-1);
    Atil = Atil_Solver([x1, x2, x3, x4]);
    
    Ftil = eye(n) + dt*Atil;
    pertX(:, ii) = Ftil*pertX(:, ii-1);
end

LinX(1, :) = nomCon(1,:) + pertX(1, :);
LinX(2, :) = nomCon(2,:) + pertX(2, :);
LinX(3, :) = nomCon(3,:) + pertX(3, :);
LinX(4, :) = nomCon(4,:) + pertX(4, :);

%% perturbations plot
figure()
subplot(4, 1, 1)
plot(Time_out, pertX(1, :))
xlabel('Time [s]')
ylabel('\deltaX [km]')

subplot(4, 1, 2)
plot(Time_out, pertX(2, :))
xlabel('Time [s]')
ylabel('\deltaXdot [km/s]')

subplot(4, 1, 3)
plot(Time_out, pertX(3, :))
xlabel('Time [s]')
ylabel('\deltaY [km]')

subplot(4, 1, 4)
plot(Time_out, pertX(4, :))
xlabel('Time [s]')
ylabel('\deltaYdot [km/s]')

sgtitle('Linearized Approx Perturbations vs Time')

% state plot
figure()
subplot(4, 1, 1)
plot(Time_out, LinX(1, :))
xlabel('Time [s]')
ylabel('X [km]')

subplot(4, 1, 2)
plot(Time_out, LinX(2, :))
xlabel('Time [s]')
ylabel('Xdot [km/s]')

subplot(4, 1, 3)
plot(Time_out, LinX(3, :))
xlabel('Time [s]')
ylabel('Y [km]')

subplot(4, 1, 4)
plot(Time_out, LinX(4, :))
xlabel('Time [s]')
ylabel('Ydot [km/s]')

sgtitle('States vs Time, Linearized Approximate Dynamics Soluiton')

%% Everything below here is wrong

% Measurement Lin
% for ii = 1:12
%     rhoLinPert(:, ii) = sqrt((pertX(1,:)' - TS_X(:, ii)).^2 + (pertX(3,:)' - TS_Y(:, ii)).^2);
%     rho_dotLinPert(:, ii) = ((pertX(1,:)' - TS_X(:, ii)).*(pertX(2,:)' - TS_Xdot(:, ii)) + (pertX(3,:)' - TS_Y(:, ii)).*(pertX(4,:)' - TS_Ydot(:, ii)))./rhoLinPert(:, ii);
%     phiLinPert(:, ii) = atan2((pertX(3,:)' - TS_Y(:, ii)), (pertX(1,:)' - TS_X(:, ii)));
% end

for j = 1:12
%     pertY(1,1) = sqrt((pertX(1,1)' - TS_X(1, j)).^2 + (pertX(3,1)' - TS_Y(1, j)).^2);
%     pertY(2,1) = ((pertX(1,1)' - TS_X(1, j)).*(pertX(2,1)' - TS_Xdot(1, j)) + (pertX(3,1)' - TS_Y(1, j)).*(pertX(4,1)' - TS_Ydot(1, j)))./pertY(1, 1);
%     pertY(3,1) = atan2((pertX(3,1)' - TS_Y(1, j)), (pertX(1,1)' - TS_X(1, j)));
    for ii = 1:1401
        x1 = nomCon(1,ii);
        x2 = nomCon(2,ii);
        x3 = nomCon(3,ii);
        x4 = nomCon(4,ii);
        z1 = TS_X(ii,1);
        z2 = TS_Xdot(ii,1);
        z3 = TS_Y(ii,1);
        z4 = TS_Ydot(ii,1);
        Cnom = subs(C);

        H = Cnom;
        pertY(:, ii) = H*pertX(:, ii);
    end
    
    rhoLinPert(:,j) = pertY(1,:)';
    rho_dotLinPert(:,j) = pertY(2,:)';
    phiLinPert(:,j) = pertY(3,:)';

    rhoLin(:, j) = sqrt((nomCon(1,:)' - TS_X(:, j)).^2 + (nomCon(3,:)' - TS_Y(:, j)).^2);
    rho_dotLin(:, j) = ((nomCon(1,:)' - TS_X(:, j)).*(nomCon(2,:)' - TS_Xdot(:, j)) + (nomCon(3,:)' - TS_Y(:, j)).*(nomCon(4,:)' - TS_Ydot(:, j)))./rhoLin(:, j);
    phiLin(:, j) = atan2((nomCon(3,:)' - TS_Y(:, j)), (nomCon(1,:)' - TS_X(:, j)));
%     rhoLin(:, j) = sqrt((LinX(1,:)' - TS_X(:, j)).^2 + (LinX(3,:)' - TS_Y(:, j)).^2);
%     rho_dotLin(:, j) = ((LinX(1,:)' - TS_X(:, j)).*(LinX(2,:)' - TS_Xdot(:, j)) + (LinX(3,:)' - TS_Y(:, j)).*(LinX(4,:)' - TS_Ydot(:, j)))./rhoLin(:, j);
%     phiLin(:, j) = atan2((LinX(3,:)' - TS_Y(:, j)), (LinX(1,:)' - TS_X(:, j)));
end

rhoLinNom = rhoLin + rhoLinPert;
rhoDotLinNom = rho_dotLin + rho_dotLinPert;
phiLinNom = phiLin + phiLinPert;

figure
hold on
tl = tiledlayout(4,1);
title(tl, "Approximate Linearized Model Data Simulation");
xlabel(tl, "Time (secs)");
nexttile
hold on
DEBUG = 0;
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
%     vis_index = find((phiLinNom(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiLinNom(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
%             (phiLinNom(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiLinNom(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
%             (phiLinNom(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiLinNom(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    if(DEBUG == 1)
%         plot(Time_out, pi/2 + thetaBound1Pos(:,ii));
% %         hold on
%         plot(Time_out, -pi/2 + thetaBound1Neg(:,ii));
        plot(Time_out, phiLin(:,ii));
        plot(Time_out, theta_TS(:,ii));
        scatter(Time_out(vis_index), phiLin(vis_index, ii)) 
        yline(pi);
        yline(-pi);
    end
     scatter(Time_out(vis_index), rhoLinNom(vis_index,ii));
    ylabel('rho^i (km)');
end

nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), rhoDotLinNom(vis_index,ii));
    ylabel('rhodot^i (km/s)');
end
nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), phiLinNom(vis_index,ii));
    ylabel('\phi^i (rads)');
end
nexttile
hold on
for ii = 1:12
    vis_index = find((phiCompare(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiCompare(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiCompare(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), visibleStation(vis_index,ii));
    ylabel('Visible Station ID');
end


%%%^Will code as off 12/6/21 6:40 pm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Provided Ydata

Gam = [0, 0; 1, 0; 0 0; 0 1];
Omega = dt*Gam;

ydata_TS_ID = NaN*ones(2, 1401);
ydata_data = NaN*ones(2*p, 1401);
c = NaN*ones(1, 1401);

for ii = 1:1401
    if ~isempty(ydata{ii})
        [~, c(ii)] = size(ydata{ii});
        if c(ii) == 1
            ydata_TS_ID(1, ii) = ydata{ii}(4, 1);
            ydata_data(1:3, ii) = ydata{ii}(1:3, 1);
        elseif c(ii) == 2
            ydata_TS_ID(1, ii) = ydata{ii}(4, 1);
            ydata_data(1:3, ii) = ydata{ii}(1:3, 1);
            ydata_TS_ID(2, ii) = ydata{ii}(4, 2);
            ydata_data(4:6, ii) = ydata{ii}(1:3, 2);
        end
    else
        ydata_TS_ID(:, ii) = NaN*ones(2, 1);
        ydata_data(:, ii) = NaN*ones(6, 1);
    end

end

ydata_TS_ID(1, 1) = 1;

TS_state(:, :, 1) = TS_X;
TS_state(:, :, 2) = TS_Xdot;
TS_state(:, :, 3) = TS_Y;
TS_state(:, :, 4) = TS_Ydot;

figure()
subplot(4, 1, 1)
hold on
for ii = 1:12
    [~, ID_index_c] = find(ydata_TS_ID(1, :) == ii);
    scatter(tvec(ID_index_c), ydata_data(1, ID_index_c), [], ColorSet(ii, :))
    [~, ID_index_c] = find(ydata_TS_ID(2, :) == ii);
    scatter(tvec(ID_index_c), ydata_data(4, ID_index_c), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('\rho [km]')

subplot(4, 1, 2)
hold on
for ii = 1:12
    [~, ID_index_c] = find(ydata_TS_ID(1, :) == ii);
    scatter(tvec(ID_index_c), ydata_data(2, ID_index_c), [], ColorSet(ii, :))
    [~, ID_index_c] = find(ydata_TS_ID(2, :) == ii);
    scatter(tvec(ID_index_c), ydata_data(5, ID_index_c), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('\rhodot [km/s]')

subplot(4, 1, 3)
hold on
for ii = 1:12
    [~, ID_index_c] = find(ydata_TS_ID(1, :) == ii);
    scatter(tvec(ID_index_c), ydata_data(3, ID_index_c), [], ColorSet(ii, :))
    [~, ID_index_c] = find(ydata_TS_ID(2, :) == ii);
    scatter(tvec(ID_index_c), ydata_data(6, ID_index_c), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('\phi [rad]')

subplot(4, 1, 4)
hold on
for ii = 1:12
    [~, ID_index_c] = find(ydata_TS_ID(1, :) == ii);
    scatter(tvec(ID_index_c), ydata_TS_ID(1, ID_index_c), [], ColorSet(ii, :))
    [~, ID_index_c] = find(ydata_TS_ID(2, :) == ii);
    scatter(tvec(ID_index_c), ydata_TS_ID(2, ID_index_c), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('Visible Station ID')

sgtitle('Provided ydata')

%% Part II - LKF Tunning
N = 100;

for jj = 1:N
    % TMT 
    p1 = unifrnd(-0.001, 0.001);
    p2 = unifrnd(-0.001, 0.001);
    p3 = unifrnd(-0.001, 0.001);
    p4 = unifrnd(-0.001, 0.001);
    perts = [p1, p2, p3, p4];
    
    MC_Initial_State = perts + initCon;

    TMT_X(1) = MC_Initial_State(1);
    TMT_Xdot(1) = MC_Initial_State(2);
    TMT_Y(1) = MC_Initial_State(3);
    TMT_Ydot(1) = MC_Initial_State(4); 
    
    v = mvnrnd([0, 0, 0], Rtrue)';
    TMT_y_NL_out(:, 1) = StatOD_NLMeasurement([TMT_X(1), TMT_Xdot(1), TMT_Y(1), TMT_Ydot(1)], [TS_X(1, 1), TS_Xdot(1, 1), TS_Y(1, 1), TS_Ydot(1, 1)]);
    TMT_y_NL_noise_out(:, 1) = TMT_y_NL_out(:, 1) + v;
    TMT_ydata(1) = {[TMT_y_NL_noise_out(:, 1); ydata_TS_ID(1)]};
    TS_ID = NaN*ones(1401, 1);
    TS_ID(1) = ydata_TS_ID(1);
    
    
    for ii = 1:1400
        % Noise
        w = mvnrnd([0, 0], Qtrue);
        v = mvnrnd([0, 0, 0], Rtrue)';
        
        % State
        ODE45_InitialState = [TMT_X(ii), TMT_Xdot(ii), TMT_Y(ii), TMT_Ydot(ii), w(1), w(2)];
        [~, TMT_test] = ode45(@(Time, State) StatODNL_noise_ODE(Time, State), [tvec(ii) tvec(ii+1)], ODE45_InitialState, options);
        
        TMT_X(ii+1) = TMT_test(end, 1);
        TMT_Xdot(ii+1) = TMT_test(end, 2);
        TMT_Y(ii+1) = TMT_test(end, 3);
        TMT_Ydot(ii+1) = TMT_test(end, 4);
        
        % Measurment
        for kk = 1:12 % compute measurements for each ground station 
            yi = StatOD_NLMeasurement([TMT_X(ii+1), TMT_Xdot(ii+1), TMT_Y(ii+1), TMT_Ydot(ii+1)], [TS_X(ii+1, kk), TS_Xdot(ii+1, kk), TS_Y(ii+1, kk), TS_Ydot(ii+1, kk)]);
            TMT_y_ALL(ii, kk, 1) = yi(1) + v(1); % rho
            TMT_y_ALL(ii, kk, 2) = yi(2) + v(2); % rhodot
            TMT_y_ALL(ii, kk, 3) = yi(3) + v(3); % phi
        end
    end
    phiCompare = TMT_y_ALL(:, :, 3);
    
    for kk = 1:12 % compute the current visible ground station   
        vis_index = find((phiCompare(:, kk) <= (pi/2 + thetaCompare(2:end, kk)) & phiCompare(:, kk) >= (-pi/2 + thetaCompare(2:end, kk))) | ...
            (phiCompare(:, kk) <= (pi/2 + thetaBound1Pos(2:end, kk)) & phiCompare(:, kk) >= (-pi/2 + thetaBound1Neg(2:end, kk))) | ...
            (phiCompare(:, kk) <= (pi/2 + thetaBound2Pos(2:end, kk)) & phiCompare(:, kk) >= (-pi/2 + thetaBound2Neg(2:end, kk))));
    
        TMT_y_NL_noise_out(1, vis_index+1) = TMT_y_ALL(vis_index, kk, 1);
        TMT_y_NL_noise_out(2, vis_index+1) = TMT_y_ALL(vis_index, kk, 2);
        TMT_y_NL_noise_out(3, vis_index+1) = TMT_y_ALL(vis_index, kk, 3);
        
        TS_ID(vis_index+1) = repmat(kk, length(vis_index), 1);
        
        
    end
    
    for ii = 2:1401
        TMT_ydata(ii) = {[TMT_y_NL_noise_out(:, ii); TS_ID(ii)]};
    end
    
    TMT_State = [TMT_X; TMT_Xdot; TMT_Y; TMT_Ydot];

    % NEES and NIS
    Q_LKF = 50*Qtrue;
    R_LKF = Rtrue;
    
    dx0 = perts;
    P0 = diag(dx0);
    
    [P, dx, x_stds, eytil, S] = LKF_StatOD(dx0, P0, TMT_ydata, dt, Q_LKF, R_LKF, Gam, TS_state, nomCon');
    
    x = dx.pos + nomCon;
    ex = TMT_State - x;
    for ii = 1:1401
        Ex(jj, ii) = ex(:, ii)'*(P.pos(:, :, ii))^-1*ex(:, ii);
        Ey(jj, ii) = eytil(1:3, ii)'*(S(1:3, 1:3, ii))^-1*eytil(1:3, ii);
    end
    
end

%%

figure()
subplot(4, 1, 1)
plot(tvec, TMT_State(1, :), 'k')
hold on
plot(tvec, nomCon(1, :), 'b')
plot(tvec, x(1, :), 'r')
hold off
xlabel('Time [s]')
ylabel('X [km]')
legend('TMT', 'Nominal', 'Estimated')


figure()
subplot(4, 1, 1)
plot(tvec, ex(1, :), 'k')
hold on
plot(tvec, 2*x_stds(1, :), 'r')
plot(tvec, -2*x_stds(1, :), 'r')
hold off
xlabel('Time [s]')
ylabel('X Error [km]')

subplot(4, 1, 2)
plot(tvec, ex(2, :), 'k')
hold on
plot(tvec, 2*x_stds(2, :), 'r')
plot(tvec, -2*x_stds(2, :), 'r')
hold off
xlabel('Time [s]')
ylabel('Xdot Error [km/s]')

subplot(4, 1, 3)
plot(tvec, ex(3, :), 'k')
hold on
plot(tvec, 2*x_stds(3, :), 'r')
plot(tvec, -2*x_stds(3, :), 'r')
hold off
xlabel('Time [s]')
ylabel('Y Error [km]')

subplot(4, 1, 4)
plot(tvec, ex(4, :), 'k')
hold on
plot(tvec, 2*x_stds(4, :), 'r')
plot(tvec, -2*x_stds(4, :), 'r')
hold off
xlabel('Time [s]')
ylabel('Ydot Error [km/s]')

% Consitancy Plots
Ex_mean = mean(Ex);
Ey_mean = mean(Ey);

alpha = 0.05;
r1 = chi2inv(alpha/2, N*n)/N;
r2 = chi2inv(1-alpha/2, N*n)/N;

figure() % NEES
scatter(tvec, Ex_mean)
hold on
plot(tvec, repmat(r1, 1401, 1), 'r--')
plot(tvec, repmat(r2, 1401, 1), 'r--')
hold off
ylim([2 6])
xlabel('Time [s]')
ylabel('Mean \epsilon_x')
title('LKF NEES Plot')
legend('NEES @ t_k', 'r_1 Bound', 'r_2 Bound')
grid on
set(gca, 'FontSize', 14) 

alpha = 0.05;
r1 = chi2inv(alpha/2, N*p)/N;
r2 = chi2inv(1-alpha/2, N*p)/N;

figure() % NIS
scatter(tvec, Ey_mean)
hold on
plot(tvec, repmat(r1, 1401, 1), 'r--')
plot(tvec, repmat(r2, 1401, 1), 'r--')
hold off
ylim([1 5])
xlabel('Time [s]')
ylabel('Mean \epsilon_y')
title('LKF NIS Plot')
legend('NIS @ t_k', 'r_1 Bound', 'r_2 Bound')
grid on
set(gca, 'FontSize', 14) 

%% Part II - EKF Tunning

N = 50;

for jj = 1:N
    % TMT 
    p1 = unifrnd(-0.25, 0.25);
    p2 = unifrnd(-0.05, 0.05);
    p3 = unifrnd(-0.25, 0.25);
    p4 = unifrnd(-0.05, 0.05);
    perts = [p1, p2, p3, p4];
    
    MC_Initial_State = perts + initCon;

    TMT_X(1) = MC_Initial_State(1);
    TMT_Xdot(1) = MC_Initial_State(2);
    TMT_Y(1) = MC_Initial_State(3);
    TMT_Ydot(1) = MC_Initial_State(4); 
    
    v = mvnrnd([0, 0, 0], Rtrue)';
    TMT_y_NL_out(:, 1) = StatOD_NLMeasurement([TMT_X(1), TMT_Xdot(1), TMT_Y(1), TMT_Ydot(1)], [TS_X(1, 1), TS_Xdot(1, 1), TS_Y(1, 1), TS_Ydot(1, 1)]);
    TMT_y_NL_noise_out(:, 1) = TMT_y_NL_out(:, 1) + v;
    TMT_ydata(1) = {[TMT_y_NL_noise_out(:, 1); ydata_TS_ID(1)]};
    TS_ID = NaN*ones(1401, 1);
    TS_ID(1) = ydata_TS_ID(1);
    
    
    for ii = 1:1400
        % Noise
        w = mvnrnd([0, 0], Qtrue);
        v = mvnrnd([0, 0, 0], Rtrue)';
        
        % State
        ODE45_InitialState = [TMT_X(ii), TMT_Xdot(ii), TMT_Y(ii), TMT_Ydot(ii), w(1), w(2)];
        [~, TMT_test] = ode45(@(Time, State) StatODNL_noise_ODE(Time, State), [tvec(ii) tvec(ii+1)], ODE45_InitialState, options);
        
        TMT_X(ii+1) = TMT_test(end, 1);
        TMT_Xdot(ii+1) = TMT_test(end, 2);
        TMT_Y(ii+1) = TMT_test(end, 3);
        TMT_Ydot(ii+1) = TMT_test(end, 4);
        
        % Measurment
        for kk = 1:12 % compute measurements for each ground station 
            yi = StatOD_NLMeasurement([TMT_X(ii+1), TMT_Xdot(ii+1), TMT_Y(ii+1), TMT_Ydot(ii+1)], [TS_X(ii+1, kk), TS_Xdot(ii+1, kk), TS_Y(ii+1, kk), TS_Ydot(ii+1, kk)]);
            TMT_y_ALL(ii, kk, 1) = yi(1) + v(1); % rho
            TMT_y_ALL(ii, kk, 2) = yi(2) + v(2); % rhodot
            TMT_y_ALL(ii, kk, 3) = yi(3) + v(3); % phi
        end
    end
    phiCompare = TMT_y_ALL(:, :, 3);
    
    for kk = 1:12 % compute the current visible ground station   
        vis_index = find((phiCompare(:, kk) <= (pi/2 + thetaCompare(2:end, kk)) & phiCompare(:, kk) >= (-pi/2 + thetaCompare(2:end, kk))) | ...
            (phiCompare(:, kk) <= (pi/2 + thetaBound1Pos(2:end, kk)) & phiCompare(:, kk) >= (-pi/2 + thetaBound1Neg(2:end, kk))) | ...
            (phiCompare(:, kk) <= (pi/2 + thetaBound2Pos(2:end, kk)) & phiCompare(:, kk) >= (-pi/2 + thetaBound2Neg(2:end, kk))));
    
        TMT_y_NL_noise_out(1, vis_index+1) = TMT_y_ALL(vis_index, kk, 1);
        TMT_y_NL_noise_out(2, vis_index+1) = TMT_y_ALL(vis_index, kk, 2);
        TMT_y_NL_noise_out(3, vis_index+1) = TMT_y_ALL(vis_index, kk, 3);
        
        TS_ID(vis_index+1) = repmat(kk, length(vis_index), 1);
        
        
    end
    
    for ii = 2:1401
        TMT_ydata(ii) = {[TMT_y_NL_noise_out(:, ii); TS_ID(ii)]};
    end
    
    TMT_State = [TMT_X; TMT_Xdot; TMT_Y; TMT_Ydot];

    % NEES and NIS
    Q_EKF = 0.95*Qtrue;
    R_EKF = Rtrue;
    
    x0 = initCon;
    P0 = 10*eye(n);
    
    [P, x, x_stds, eytil, S] = EKF_StatOD(x0, P0, TMT_ydata, dt, tvec, Q_EKF, R_EKF, Gam, TS_state);
    
    ex = TMT_State - x.pos;
    for ii = 1:1401
        Ex(jj, ii) = ex(:, ii)'*(P.pos(:, :, ii))^-1*ex(:, ii);
        Ey(jj, ii) = eytil(1:3, ii)'*(S(1:3, 1:3, ii))^-1*eytil(1:3, ii);
    end
    
end

figure() % state estimation errors
subplot(4, 1, 1)
plot(tvec, ex(1, :), 'k')
hold on
plot(tvec, 2*x_stds(1, :), 'r--')
plot(tvec, -2*x_stds(1, :), 'r--')
hold off
ylim([-0.5 0.5])
xlabel('Time [s]')
ylabel('X Error [km]')
legend('State Estimation Error', '2\sigma Bounds')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 2)
plot(tvec, ex(2, :), 'k')
hold on
plot(tvec, 2*x_stds(2, :), 'r--')
plot(tvec, -2*x_stds(2, :), 'r--')
hold off
ylim([-0.005 0.005])
xlabel('Time [s]')
ylabel('Xdot Error [km/s]')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 3)
plot(tvec, ex(3, :), 'k')
hold on
plot(tvec, 2*x_stds(3, :), 'r--')
plot(tvec, -2*x_stds(3, :), 'r--')
hold off
ylim([-0.6 0.6])
xlabel('Time [s]')
ylabel('Y Error [km]')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 4)
plot(tvec, ex(4, :), 'k')
hold on
plot(tvec, 2*x_stds(4, :), 'r--')
plot(tvec, -2*x_stds(4, :), 'r--')
hold off
ylim([-0.005 0.005])
xlabel('Time [s]')
ylabel('Ydot Error [km/s]')
grid on
set(gca, 'FontSize', 14)

sgtitle('States Estimation Error vs Time - EKF')

figure() % innovations
subplot(3, 1, 1)
scatter(tvec, eytil(1, :))
xlabel('Time [s]')
ylabel('\rho Error [km]');
grid on
set(gca, 'FontSize', 14)
ylim([-0.5 0.5])

subplot(3, 1, 2)
scatter(tvec, eytil(2, :))
xlabel('Time [s]')
ylabel('\rhodot Error [km/s]');
grid on
set(gca, 'FontSize', 14)

subplot(3, 1, 3)
scatter(tvec, eytil(3, :))
xlabel('Time [s]')
ylabel('\phi Error [rad]');
grid on
set(gca, 'FontSize', 14)

sgtitle('Inovations vs Time - EKF')

% NEES and NIS Plots

Ex_mean = mean(Ex);
Ey_mean = mean(Ey);

alpha = 0.05;
r1 = chi2inv(alpha/2, N*n)/N;
r2 = chi2inv(1-alpha/2, N*n)/N;

figure() % NEES
scatter(tvec, Ex_mean)
hold on
plot(tvec, repmat(r1, 1401, 1), 'r--')
plot(tvec, repmat(r2, 1401, 1), 'r--')
hold off
ylim([2 6])
xlabel('Time [s]')
ylabel('Mean \epsilon_x')
title('EKF NEES Plot')
legend('NEES @ t_k', 'r_1 Bound', 'r_2 Bound')
grid on
set(gca, 'FontSize', 14)

r1 = chi2inv(alpha/2, N*p)/N;
r2 = chi2inv(1-alpha/2, N*p)/N;

figure() % NIS
scatter(tvec, Ey_mean);
hold on
plot(tvec, repmat(r1, 1401, 1), 'r--')
plot(tvec, repmat(r2, 1401, 1), 'r--')
hold off
ylim([1 5])
xlabel('Time [s]')
ylabel('Mean \epsilon_y')
title('EKF NIS Plot')
legend('NIS @ t_k', 'r_1 Bound', 'r_2 Bound')
grid on
set(gca, 'FontSize', 14)

figure() % TMT state plots
subplot(4, 1, 1)
plot(tvec, TMT_X, 'k')
hold on
plot(tvec, nomCon(1, :), 'r')
hold off
xlabel('Time [s]')
ylabel('TMT X [km]')
legend('TMT State', 'Nominal Trajectory')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 2)
plot(tvec, TMT_Xdot, 'k')
hold on
plot(tvec, nomCon(2, :), 'r')
hold off
xlabel('Time [s]')
ylabel('TMT Xdot [km/s]')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 3)
plot(tvec, TMT_Y, 'k')
hold on
plot(tvec, nomCon(3, :), 'r')
hold off
xlabel('Time [s]')
ylabel('TMT Y [km]')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 4)
plot(tvec, TMT_Ydot, 'k')
hold on
plot(tvec, nomCon(4, :), 'r')
hold off
xlabel('Time [s]')
ylabel('TMT Ydot [km/s]')
grid on
set(gca, 'FontSize', 14)

sgtitle('TMT Simulated States vs Time')

figure() % TMT Measurement Plots
subplot(4, 1, 1)
hold on
for ii = 1:12
    ID_index = find(TS_ID == ii);
    scatter(tvec(ID_index), TMT_y_NL_noise_out(1, ID_index), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('TMT \rho [km]')
grid on 
set(gca, 'FontSize', 14)

subplot(4, 1, 2)
hold on
for ii = 1:12
    ID_index = find(TS_ID == ii);
    scatter(tvec(ID_index), TMT_y_NL_noise_out(2, ID_index), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('TMT \rhodot [km/s]')
grid on 
set(gca, 'FontSize', 14)

subplot(4, 1, 3)
hold on
for ii = 1:12
    ID_index = find(TS_ID == ii);
    scatter(tvec(ID_index), TMT_y_NL_noise_out(3, ID_index), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('TMT \phi [rad]')
grid on 
set(gca, 'FontSize', 14) 

subplot(4, 1, 4)
hold on
for ii = 1:12
    ID_index = find(TS_ID == ii);
    scatter(tvec(ID_index), TS_ID(ID_index), [], ColorSet(ii, :))
end
hold off
xlabel('Time [s]')
ylabel('Visible Station ID')
grid on 
set(gca, 'FontSize', 14)

sgtitle('TMT Simulated Measurements vs Time')

%% EKF Implementation

x0 = Initial_States;
P0 = 10*eye(n);

[~, x, x_stds, ~, ~] = EKF_StatOD(x0, P0, ydata, dt, tvec, Q_EKF, R_EKF, Gam, TS_state);

figure()
subplot(4, 1, 1)
plot(tvec, x.pos(1, :), 'k')
hold on
plot(tvec, x.pos(1, :) + 2*x_stds(1, :), 'r--')
plot(tvec, x.pos(1, :) - 2*x_stds(1, :), 'r--')
hold off
xlabel('Time [s]')
ylabel('Estimated X [km]')
legend('Estimated State', '2\sigma')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 2)
plot(tvec, x.pos(2, :), 'k')
hold on
plot(tvec, x.pos(2, :) + 2*x_stds(2, :), 'r--')
plot(tvec, x.pos(2, :) - 2*x_stds(2, :), 'r--')
hold off
xlabel('Time [s]')
ylabel('Estimated Xdot [km/s]')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 3)
plot(tvec, x.pos(3, :), 'k')
hold on
plot(tvec, x.pos(3, :) + 2*x_stds(3, :), 'r--')
plot(tvec, x.pos(3, :) - 2*x_stds(3, :), 'r--')
hold off
xlabel('Time [s]')
ylabel('Estimated Y [km]')
grid on
set(gca, 'FontSize', 14)

subplot(4, 1, 4)
plot(tvec, x.pos(4, :), 'k')
hold on
plot(tvec, x.pos(4, :) + 2*x_stds(4, :), 'r--')
plot(tvec, x.pos(4, :) - 2*x_stds(4, :), 'r--')
hold off
xlabel('Time [s]')
ylabel('Estimated Ydot [km/s]')
grid on
set(gca, 'FontSize', 14)

sgtitle('Implemented EKF Estimated States vs Time')

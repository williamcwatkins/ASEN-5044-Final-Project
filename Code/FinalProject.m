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
v = sqrt(u/r0); % orbital velocity
circ = 2*pi*r0; % circumference of the orbit
T = circ / v; % Orbital period

n = 4; % number of states
m = 2; % number of inputs, also happens to be number of distrubances
p = 3; % number of measurements
j = 6; % number of measurement stations

dt = 10; % step size [s]
tspan = 0:dt:14000;
initCon = [r0, 0, 0, r0*sqrt(u/r0^3)]; % initial state point - LTV sys, so
% have to linearize about a nominal trajectory, see below

DEBUG = 0; % Debug boolean

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
for i = 1:length(tspan)
    nomCon(:,i) = [r0 * cos(sqrt(u/r0^3)*tspan(i)); -r0 * sin(sqrt(u/r0^3)*tspan(i))*sqrt(u/r0^3);
            r0*sin(sqrt(u/r0^3)*tspan(i)); r0*cos(sqrt(u/r0^3)*tspan(i))*sqrt(u/r0^3)];
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

% Visible Tracking Stations %%%%%%%%%%%%%% there is a wrapping issue
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
    Atil = subs(A);
    
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

% measurement and vis plot

%% Everything below here is wrong

% Measurement Lin
% for ii = 1:12
%     rhoLinPert(:, ii) = sqrt((pertX(1,:)' - TS_X(:, ii)).^2 + (pertX(3,:)' - TS_Y(:, ii)).^2);
%     rho_dotLinPert(:, ii) = ((pertX(1,:)' - TS_X(:, ii)).*(pertX(2,:)' - TS_Xdot(:, ii)) + (pertX(3,:)' - TS_Y(:, ii)).*(pertX(4,:)' - TS_Ydot(:, ii)))./rhoLinPert(:, ii);
%     phiLinPert(:, ii) = atan2((pertX(3,:)' - TS_Y(:, ii)), (pertX(1,:)' - TS_X(:, ii)));
% end

for ii = 1:12
    rhoLinPert(:, ii) = sqrt((pertX(1,:)').^2 + (pertX(3,:)').^2);
    rho_dotLinPert(:, ii) = ((pertX(1,:)').*(pertX(2,:)') + (pertX(3,:)').*(pertX(4,:)'))./rhoLinPert(:, ii);
    phiLinPert(:, ii) = atan2((pertX(3,:)'), (pertX(1,:)'));
end

for ii = 1:12
    rhoLin(:, ii) = sqrt((LinX(1,:)' - TS_X(:, ii)).^2 + (LinX(3,:)' - TS_Y(:, ii)).^2);
    rho_dotLin(:, ii) = ((LinX(1,:)' - TS_X(:, ii)).*(LinX(2,:)' - TS_Xdot(:, ii)) + (LinX(3,:)' - TS_Y(:, ii)).*(LinX(4,:)' - TS_Ydot(:, ii)))./rhoLin(:, ii);
    phiLin(:, ii) = atan2((LinX(3,:)' - TS_Y(:, ii)), (LinX(1,:)' - TS_X(:, ii)));
end

rhoLinNom(:,:) = rho + rhoLinPert;
rhoDotLinNom = rho_dotLinPert + rho_dot;
phiLinNom = phi + phiLinPert;

figure
hold on
tl = tiledlayout(4,1);
title(tl, "Approximate Linearized Model Data Simulation");
xlabel(tl, "Time (secs)");
nexttile
hold on
DEBUG = 0;
for ii = 1:12
    vis_index = find((phiLinNom(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiLinNom(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiLinNom(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiLinNom(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiLinNom(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiLinNom(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
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
     scatter(Time_out(vis_index), rhoLin(vis_index,ii));
    ylabel('rho^i (km)');
end

nexttile
hold on
for ii = 1:12
    vis_index = find((phiLin(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiLin(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiLin(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), rho_dotLin(vis_index,ii));
    ylabel('rhodot^i (km/s)');
end
nexttile
hold on
for ii = 1:12
    vis_index = find((phiLin(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiLin(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiLin(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), phiLinNom(vis_index,ii));
    ylabel('\phi^i (rads)');
end
nexttile
hold on
for ii = 1:12
    vis_index = find((phiLin(:, ii) <= (pi/2 + thetaCompare(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaCompare(:, ii))) | ...
            (phiLin(:, ii) <= (pi/2 + thetaBound1Pos(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaBound1Neg(:, ii))) | ...
            (phiLin(:, ii) <= (pi/2 + thetaBound2Pos(:, ii)) & phiLin(:, ii) >= (-pi/2 + thetaBound2Neg(:, ii))));
    scatter(Time_out(vis_index), visibleStation(vis_index,ii));
    ylabel('Visible Station ID');
end

%% Part II - EKF

 










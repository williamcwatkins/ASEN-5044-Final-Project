%% Admin
%{

Names: Corey LePine and William Watkins 
Professor: McMahon
Class: ASEN 5044 Stat Est for Dyn Sys
Date: December 14, 2021
Final Project: Statistical Orbit Determination

%} 
%% Housekeeping
clc; clear; close;

%% Constants
load('orbitdeterm_finalproj_KFdata.mat')
u = 398600; % Earth's standard gravitational paremters [km^3/s^2]
r0 = 6678; % Nominal orbit radius [km]
Re = 6378; % Uniform radius of Earth [km]
we = 2*pi/86400; % Consant rotation rate of Earth [rad/s]

n = 4; % number of states
m = 2; % number of inputs, also happens to be number of distrubances
p = 3; % number of measurements

dt = 10; % step size [s]
NomCon = [r0, 0, 0, r0*sqrt(u/r0^3)]; % nominal state point

%% Part 1

% 1a) Find the CT Jacobian Matracies
syms x1 x2 x3 x4 mu u1 u2 w1 w2 z1 z2 z3 z4 % zi is the states of the ground stations

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
mu = u;
x1 = NomCon(1); 
x2 = NomCon(2);
x3 = NomCon(3);
x4 = NomCon(4);
Atil = subs(A);
    
%% 1c) Nonlinear Dynamics
% State NL 
Rel_Tol = 1e-13;
Abs_Tol = Rel_Tol;
options = odeset('Stats', 'off', 'RelTol', Rel_Tol, 'AbsTol', Abs_Tol);

tspan = 0:dt:14000;
perts = [0, 0.075, 0, -0.021];
Initial_States = perts + NomCon;

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

% Measurement NL
for ii = 1:12
    rho(:, ii) = sqrt((State_X - TS_X(:, ii)).^2 + (State_Y - TS_Y(:, ii)).^2);
    rho_dot(:, ii) = ((State_X - TS_X(:, ii)).*(State_Xdot - TS_Xdot(:, ii)) + (State_Y - TS_Y(:, ii)).*(State_Ydot - TS_Ydot(:, ii)))./rho(:, ii);
    phi(:, ii) = atan2((State_Y - TS_Y(:, ii)), (State_X - TS_X(:, ii)));
end

% Visibile Tracking Stations %%%%%%%%%%%%%% there is a wrapping issue
figure
for ii = 1:12
    vis_index = find(phi(:, ii) <= (pi/2 + theta_TS(:, ii)) & phi(:, ii) >= (-pi/2 + theta_TS(:, ii)));
    hold on
    scatter(Time_out(vis_index), rho(vis_index, ii)) 
end

%% 1c) Linearized DT
pertX(:, 1) = perts';
for ii = 2:1401
    x1 = State_X(ii-1);
    x2 = State_Xdot(ii-1);
    x3 = State_Y(ii-1);
    x4 = State_Ydot(ii-1);
    
    Atil = subs(A);
    
    Ftil = eye(n) + dt*Atil;
    pertX(:, ii) = Ftil*pertX(:, ii-1);
end

LinX(1, :) = State_X' + pertX(1, :);
LinX(2, :) = State_Xdot' + pertX(2, :);
LinX(3, :) = State_Y' + pertX(3, :);
LinX(4, :) = State_Ydot' + pertX(4, :);

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


%% Part II - EKF

 










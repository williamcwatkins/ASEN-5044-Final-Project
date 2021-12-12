function [State_Derivatives] = StatODNL_noise_ODE(Time, State)

u = 398600; % Earth's standard gravitational paremters [km^3/s^2]

X = State(1);
Xdot = State(2);
Y = State(3);
Ydot = State(4);

r = sqrt(X^2 + Y^2);

w1 = State(5);
w2 = State(6);

Xddot = -u*X/r^3 + w1;
Yddot = -u*Y/r^3 + w2;

State_Derivatives = [Xdot; Xddot; Ydot; Yddot; 0; 0];

end


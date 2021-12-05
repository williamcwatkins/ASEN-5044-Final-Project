function [State_Derivatives] = StatODNL_ODE(Time, State)

u = 398600; % Earth's standard gravitational paremters [km^3/s^2]

X = State(1);
Xdot = State(2);
Y = State(3);
Ydot = State(4);

r = sqrt(X^2 + Y^2);

Xddot = -u*X/r^3;
Yddot = -u*Y/r^3;

State_Derivatives = [Xdot; Xddot; Ydot; Yddot];

end


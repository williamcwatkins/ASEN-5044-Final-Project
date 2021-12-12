function [y] = StatOD_NLMeasurement(Sat_state, TS_state)

X = Sat_state(1);
Xdot = Sat_state(2);
Y = Sat_state(3);
Ydot = Sat_state(4);

Xi = TS_state(1);
Xidot = TS_state(2);
Yi = TS_state(3);
Yidot = TS_state(4);

rho = sqrt((X - Xi)^2 + (Y - Yi)^2);
rhodot = ((X - Xi)*(Xdot - Xidot) + (Y - Yi)*(Ydot - Yidot))/sqrt((X - Xi)^2 + (Y - Yi)^2);
phi = atan2(Y - Yi, X - Xi);

y = [rho; rhodot; phi];

end


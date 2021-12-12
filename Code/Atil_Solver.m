function [Atil] = Atil_Solver(State)

mu = 398600; % Earth's standard gravitational paremters [km^3/s^2]

x1 = State(1);
% x2 = State(2);
x3 = State(3);
% x4 = State(4);

Atil = [0, 1, 0, 0;
    (3*mu*x1^2)/(x1^2 + x3^2)^(5/2) - mu/(x1^2 + x3^2)^(3/2), 0, (3*mu*x1*x3)/(x1^2 + x3^2)^(5/2), 0;
    0, 0, 0, 1;
    (3*mu*x1*x3)/(x1^2 + x3^2)^(5/2), 0, (3*mu*x3^2)/(x1^2 + x3^2)^(5/2) - mu/(x1^2 + x3^2)^(3/2), 0];

end


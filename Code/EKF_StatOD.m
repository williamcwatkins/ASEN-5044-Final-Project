function [P, x, x_stds, eytil, S] = EKF_StatOD(x0, P0, ydata, dt, tvec, Q, R, Gamma, TS_state)


% DEAL WITH MULT Measurement values!!!!!!!

% Matrix sizes and Steps
n = length(x0); % number of states
p = length(ydata{1}) - 1; % number of measurments, subtract one since it has GND station ID 
steps = length(ydata); % number of steps for problem; step 1 is the zero time vec

%mu = 398600; % Earth's standard gravitational paremters [km^3/s^2]
Omega = dt*Gamma; % Since CT Gamma matrix is LTI we can compute Omega outside EKF loop

% ODE Tolerances
Rel_Tol = 1e-13;
Abs_Tol = Rel_Tol;
options = odeset('Stats', 'off', 'RelTol', Rel_Tol, 'AbsTol', Abs_Tol);

% initilaize variables for speed
x.neg(1:n, 1) = NaN*ones(n, 1);
x.pos(:, 1) = x0;

P.neg(1:n, 1:n, 1) = NaN*ones(n);
P.pos(:, :, 1) = P0;

x_stds(1:n, 1) = sqrt(diag(P0));

eytil = NaN*ones(2*p, 1);

S = NaN*ones(2*p, 2*p, steps);

y = NaN*ones(2*p, steps);
TS_ID = NaN*ones(2, steps);
c = NaN*ones(1, steps);

% Parse ydata into own data vector and GND station IDS
for ii = 1:steps
    if ~isempty(ydata{ii})
        [~, c(ii)] = size(ydata{ii});
        if c(ii) == 1
            y(1:3, ii) = ydata{ii}(1:3);
            TS_ID(1, ii) = ydata{ii}(4);
        elseif c(ii) == 2
            y(1:3, ii) = ydata{ii}(1:3, 1);
            TS_ID(1, ii) = ydata{ii}(4, 1);
            y(4:6, ii) = ydata{ii}(1:3, 2);
            TS_ID(2, ii) = ydata{ii}(4, 2);
        end
    end
end

No_Meas_index = isnan(TS_ID(1, :));

% EKF Loop
for ii = 2:steps
   % Prediction Step
   tspan = [tvec(ii-1) tvec(ii)];
   [~, NL_state] = ode45(@(Time, State) StatODNL_ODE(Time, State), tspan, x.pos(:, ii-1)', options);
   x.neg(:, ii) = NL_state(end, :)';

   % NL and Jacobian Computaion (Part of Prediction Step)
   Atil = Atil_Solver(x.pos(:, ii-1));
   
   Ftil = eye(n) + dt*Atil;
   P.neg(:, :, ii) = Ftil*P.pos(:, :, ii-1)*Ftil' + Omega*Q*Omega';
   
   % Correction Step
   if No_Meas_index(ii) == 0 % There are Measurments
       if c(ii) == 1
           z1 = TS_state(ii, TS_ID(1, ii), 1);
           z2 = TS_state(ii, TS_ID(1, ii), 2);
           z3 = TS_state(ii, TS_ID(1, ii), 3);
           z4 = TS_state(ii, TS_ID(1, ii), 4);
           TS_stateK = [z1; z2; z3; z4];

           y_neg = StatOD_NLMeasurement(x.neg(:, ii), TS_stateK);

           Htil = Ctil_Solver(x.neg(:, ii), TS_stateK);

           eytil(1:3, ii) = y(1:3, ii) - y_neg;

           Ktil = P.neg(:, :, ii)*Htil'*(Htil*P.neg(:, :, ii)*Htil' + R)^-1;

           x.pos(:, ii) = x.neg(:, ii) + Ktil*eytil(1:3, ii);
           P.pos(:, :, ii) = (eye(n) - Ktil*Htil)*P.neg(:, :, ii);

           S(1:p, 1:p, ii) = Htil*P.neg(:, :, ii)*Htil' + R;
       
       elseif c(ii) == 2
           % first measurment
           z1 = TS_state(ii, TS_ID(1, ii), 1);
           z2 = TS_state(ii, TS_ID(1, ii), 2);
           z3 = TS_state(ii, TS_ID(1, ii), 3);
           z4 = TS_state(ii, TS_ID(1, ii), 4);
           TS_stateK = [z1; z2; z3; z4];

           y_neg_1 = StatOD_NLMeasurement(x.neg(:, ii), TS_stateK);

           Htil_1 = Ctil_Solver(x.neg(:, ii), TS_stateK);

           eytil(1:3, ii) = y(1:3, ii) - y_neg_1;
           Ktil_1 = P.neg(:, :, ii)*Htil_1'*(Htil_1*P.neg(:, :, ii)*Htil_1' + R)^-1;
           
           % second measurment
           z1 = TS_state(ii, TS_ID(2, ii), 1);
           z2 = TS_state(ii, TS_ID(2, ii), 2);
           z3 = TS_state(ii, TS_ID(2, ii), 3);
           z4 = TS_state(ii, TS_ID(2, ii), 4);
           TS_stateK = [z1; z2; z3; z4];

           y_neg_2 = StatOD_NLMeasurement(x.neg(:, ii), TS_stateK);

           Htil_2 = Ctil_Solver(x.neg(:, ii), TS_stateK);

           eytil(4:6, ii) = y(4:6, ii) - y_neg_2;
           
           Ktil_2 = P.neg(:, :, ii)*Htil_2'*(Htil_2*P.neg(:, :, ii)*Htil_2' + R)^-1;
           
           % Combined
           Pblock = blkdiag(P.neg(:, :, ii), P.neg(:, :, ii));
           Hblock = blkdiag(Htil_1, Htil_2);
           Rblock = blkdiag(R, R);
           
           x.pos(:, ii) = x.neg(:, ii) + Ktil_1*eytil(1:3, ii) + Ktil_2*eytil(4:6, ii);
           P.pos(:, :, ii) =  P.neg(:, :, ii) - Ktil_1*Htil_1*P.neg(:, :, ii) - Ktil_2*Htil_2*P.neg(:, :, ii);
           
           S(:, :, ii) = Hblock*Pblock*Hblock' + Rblock;
       end      
   else % No Measurements
       x.pos(:, ii) = x.neg(:, ii);
       P.pos(:, :, ii) = P.neg(:, :, ii);
       eytil(:, ii) = NaN*ones(2*p, 1);
   end
   
   % Standard Deviations
   x_stds(:, ii) = sqrt(diag(P.pos(:, :, ii)));
   
end

end


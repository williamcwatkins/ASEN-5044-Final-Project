function [P, dx, x_stds, eytil, S] = LKF_StatOD(dx0, P0, ydata, dt, Q, R, Gamma, TS_state, Nom_State)

% Matrix sizes and Steps
n = length(dx0); % number of states
p = length(ydata{1}) - 1; % number of measurments, subtract one since it has GND station ID 
steps = length(ydata); % number of steps for problem; step 1 is the zero time vec

Omega = dt*Gamma; % Since CT Gamma matrix is LTI we can compute Omega outside EKF loop

% initilaize variables for speed
dx.neg(1:n, 1) = NaN*ones(n, 1);
dx.pos(:, 1) = dx0;

P.neg(1:n, 1:n, 1) = NaN*ones(n);
P.pos(:, :, 1) = P0;

x_stds(1:n, 1) = sqrt(diag(P0));

eytil(1:p, 1) = NaN*ones(p, 1);

S = NaN*ones(p, p, steps);

y = NaN*ones(p, steps);
TS_ID = NaN*ones(1, steps);

% Parse ydata into own data vector and GND station IDS
for ii = 1:steps
    if ~isempty(ydata{ii})
        y(:, ii) = ydata{ii}(1:3);
        TS_ID(ii) = ydata{ii}(4);
    else
        y(:, ii) = NaN*ones(p, 1);
        TS_ID(ii) = NaN;
    end
end

% LKF Loop
for ii = 2:steps
    % Prediction Step
    Atil = Atil_Solver(Nom_State(ii-1, :));
    Ftil = eye(n) + dt*Atil;
    
    dx.neg(:, ii) = Ftil*dx.pos(:, ii-1);
    P.neg(:, :, ii) = Ftil*P.pos(:, :, ii-1)*Ftil' + Omega*Q*Omega';
    
    % Correction Step %%% Deal with multiple measurments
    if ~isnan(TS_ID(ii))
        z1 = TS_state(ii, TS_ID(ii), 1);
        z2 = TS_state(ii, TS_ID(ii), 2);
        z3 = TS_state(ii, TS_ID(ii), 3);
        z4 = TS_state(ii, TS_ID(ii), 4);
        TS_stateK = [z1; z2; z3; z4];
        
        Htil = Ctil_Solver(Nom_State(ii, :), TS_stateK);
        
        Ktil = P.neg(:, :, ii)*Htil'*(Htil*P.neg(:, :, ii)*Htil' + R)^-1;
        
        y_nom = StatOD_NLMeasurement(Nom_State(ii, :), TS_stateK);
        
        dy = y(:, ii) - y_nom;
        
        dx.pos(:, ii) = dx.neg(:, ii) + Ktil*(dy - Htil*dx.neg(:, ii));
        
        P.pos(:, :, ii) = (eye(n) - Ktil*Htil)*P.neg(:, :, ii);
        
        S(:, :, ii) = Htil*P.neg(:, :, ii)*Htil' + R;
        eytil(:, ii) = dy - Htil*dx.neg(:, ii);
       
    else
       dx.pos(:, ii) = dx.neg(:, ii);
       P.pos(:, :, ii) = P.neg(:, :, ii);
       eytil(:, ii) = NaN*ones(p, 1);
       S(:, :, ii) = NaN*ones(p);
    end
    
    % Standard Deviations
    x_stds(:, ii) = sqrt(diag(P.pos(:, :, ii))); 
    
end

end


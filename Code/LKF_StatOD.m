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

% LKF Loop
for ii = 2:steps
    % Prediction Step
    Atil = Atil_Solver(Nom_State(ii-1, :));
    Ftil = eye(n) + dt*Atil;
    
    dx.neg(:, ii) = Ftil*dx.pos(:, ii-1);
    
    P.neg(:, :, ii) = Ftil*P.pos(:, :, ii-1)*Ftil' + Omega*Q*Omega';
    
    % Correction Step %%% Deal with multiple measurments
    if No_Meas_index(ii) == 0 % There are Measurments
        if c(ii) == 1
        z1 = TS_state(ii, TS_ID(1, ii), 1);
        z2 = TS_state(ii, TS_ID(1, ii), 2);
        z3 = TS_state(ii, TS_ID(1, ii), 3);
        z4 = TS_state(ii, TS_ID(1, ii), 4);
        TS_stateK = [z1; z2; z3; z4];
        
        Htil = Ctil_Solver(Nom_State(ii, :), TS_stateK);
        
        Ktil = P.neg(:, :, ii)*Htil'*(Htil*P.neg(:, :, ii)*Htil' + R)^-1;
        
        y_nom = StatOD_NLMeasurement(Nom_State(ii, :), TS_stateK);
        
        dy = y(1:3, ii) - y_nom;
        
        dx.pos(:, ii) = dx.neg(:, ii) + Ktil*(dy - Htil*dx.neg(:, ii));
        
        P.pos(:, :, ii) = (eye(n) - Ktil*Htil)*P.neg(:, :, ii);
        
        S(1:3, 1:3, ii) = Htil*P.neg(:, :, ii)*Htil' + R;
        eytil(1:3, ii) = dy - Htil*dx.neg(:, ii);
        
        elseif c(ii) == 2
            % First measurment
            z1 = TS_state(ii, TS_ID(1, ii), 1);
            z2 = TS_state(ii, TS_ID(1, ii), 2);
            z3 = TS_state(ii, TS_ID(1, ii), 3);
            z4 = TS_state(ii, TS_ID(1, ii), 4);
            TS_stateK = [z1; z2; z3; z4];

            Htil_1 = Ctil_Solver(Nom_State(ii, :), TS_stateK);

            Ktil_1 = P.neg(:, :, ii)*Htil_1'*(Htil_1*P.neg(:, :, ii)*Htil_1' + R)^-1;

            y_nom_1 = StatOD_NLMeasurement(Nom_State(ii, :), TS_stateK);

            dy_1 = y(1:3, ii) - y_nom_1;
            
            eytil(1:3, ii) = dy_1 - Htil_1*dx.neg(:, ii);
            
            % Second measurment
            z1 = TS_state(ii, TS_ID(2, ii), 1);
            z2 = TS_state(ii, TS_ID(2, ii), 2);
            z3 = TS_state(ii, TS_ID(2, ii), 3);
            z4 = TS_state(ii, TS_ID(2, ii), 4);
            TS_stateK = [z1; z2; z3; z4];

            Htil_2 = Ctil_Solver(Nom_State(ii, :), TS_stateK);

            Ktil_2 = P.neg(:, :, ii)*Htil_2'*(Htil_2*P.neg(:, :, ii)*Htil_2' + R)^-1;

            y_nom_2 = StatOD_NLMeasurement(Nom_State(ii, :), TS_stateK);

            dy_2 = y(4:6, ii) - y_nom_2;
            
            eytil(4:6, ii) = dy_2 - Htil_2*dx.neg(:, ii);
            
            % Combined 
            Pblock = blkdiag(P.neg(:, :, ii), P.neg(:, :, ii));
            Hblock = blkdiag(Htil_1, Htil_2);
            Rblock = blkdiag(R, R);
            
            dx.pos(:, ii) = dx.neg(:, ii) + Ktil_1*(dy_1 - Htil_1*dx.neg(:, ii)) + Ktil_2*(dy_2 - Htil_2*dx.neg(:, ii));
        
            P.pos(:, :, ii) = P.neg(:, :, ii) - Ktil_1*Htil_1*P.neg(:, :, ii) - Ktil_2*Htil_2*P.neg(:, :, ii);
            
            S(:, :, ii) = Hblock*Pblock*Hblock' + Rblock;
            
        end      
    else
       dx.pos(:, ii) = dx.neg(:, ii);
       P.pos(:, :, ii) = P.neg(:, :, ii);
       eytil(:, ii) = NaN*ones(2*p, 1);
    end
    
    % Standard Deviations
    x_stds(:, ii) = sqrt(diag(P.pos(:, :, ii))); 
    
end

end


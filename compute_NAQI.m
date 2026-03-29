% compute the quantum parameter of NAQI
function [N_max,best_x] = compute_NAQI(rho_AB)

    ga_opts = optimoptions('ga', ...
    'Display', 'off', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 50, ...
    'ConstraintTolerance', 1e-8);  
    
    options = optimoptions('fmincon',...
    'Display','off',...
    'MaxFunctionEvaluations',5000,...
    'ConstraintTolerance', 1e-8,...
    'Algorithm','sqp');
    
    n = 8;
    lb = [0, 0, 0, 0, 0, 0, 0, 0];
    ub = [pi, 2*pi, pi, 2*pi, pi, 2*pi, pi, 2*pi];
    num = 50;
    
    N_max = -inf;
    for i = 1:num
        [x0, ~] = ga(@(x) objective(x, rho_AB), n, [], [], [], [], lb, ub, [], ga_opts);
        [x_opt, fval] = fmincon(@(x) -objective(x, rho_AB), x0, [], [], [], [], lb, ub,  [], options);
        candidate = -fval;
        if candidate > N_max
            N_max = candidate;
            best_x = x_opt; % 更新最优角度
%             disp(['The result of the ', num2str(i), '-th optimaization is',num2str(N_max)])
        end
    end
end





%% objective function
function N = objective(x, rho)
    H = [1;0]; V = [0;1];
    
    % Measurement basis for Bob
    M1 = {cos(x(1)/2)*H + exp(1i*x(2))*sin(x(1)/2)*V, 
          sin(x(1)/2)*H - exp(1i*x(2))*cos(x(1)/2)*V};
    M2 = {(M1{1} + M1{2})/sqrt(2), (M1{1} - M1{2})/sqrt(2)};
    M3 = {(M1{1} + 1i*M1{2})/sqrt(2), (M1{1} - 1i*M1{2})/sqrt(2)};
    

    V1 = [M1{1}, M1{2}];
    V2 = [M2{1}, M2{2}];
    V3 = [M3{1}, M3{2}];
    bases = {V1, V2, V3}; 
    
    % Measurement setting for Alice
    measurements = cell(1,3);
    for i = 1:3
        theta = x(2*i+1:2*i+2);
        op1 = cos(theta(1)/2)*H + exp(1i*theta(2))*sin(theta(1)/2)*V;
        op2 = sin(theta(1)/2)*H - exp(1i*theta(2))*cos(theta(1)/2)*V;
        measurements{i} = {op1, op2};
    end
    
    N = 0;
    for i = 1:3 
        U = bases{i}; 
        for a = 1:2  
            [prob, rho_cond] = probability_B(rho, measurements{i}{a});
            if prob > 1e-10  % ignoring events with zero probaility
                rho_trans = transform_B(rho_cond, U);
                rho_ = reduced_matrix_B(rho_trans);
                fai_val = fai_calculate(rho_);
                N = N + prob * fai_val;
            end
        end
    end
end

%%
function [prob, rho_cond] = probability_B(rho, op)

    MA = op * op'; 
    MA_full = kron(MA, eye(2)); 
    rho_cond = MA_full * rho * MA_full';
    prob = real(trace(rho_cond));  
    
    if prob > 1e-10
        rho_cond = rho_cond / prob;
    else
        rho_cond = zeros(4);
    end
end

%% 
function rho_trans = transform_B(rho, U) 
    a = [U(1,1);U(2,1)];
    b = [U(1,2);U(2,2)];
    U1 = [a,b]';
    U_full = kron(U1,U1);
    rho_trans = U_full * rho * U_full';
end

%% reduced Matrix for Alice
function rho_A = reduced_matrix_A(rho)
    % 计算B系统的约化密度矩阵
    rho_A = zeros(2);
    for i = 1:2
        idx_i = [2*i-1, 2*i];
        for j = 1:2
            idx_j = [2*j-1, 2*j];
            rho_A(i,j) = sum(diag(rho(idx_i, idx_j)));
        end
    end
    rho_A = rho_A / trace(rho_A);  % 确保归一化
end

%% reduced Matrix for Bob
function rho_B = reduced_matrix_B(rho)
    % 计算B系统的约化密度矩阵
    rho_B = zeros(2);
    for i = 1:2
        for j = 1:2
            idx_i = [i, i+2];
            idx_j = [j, j+2];
            rho_B(i,j) = sum(diag(rho(idx_i, idx_j)));
        end
    end
    rho_B = rho_B / trace(rho_B);  % 确保归一化
end

%% quantum imaginarity based on l1-norm
function f = fai_calculate(rho)
    f = sum(sum(abs(imag(rho))));
end

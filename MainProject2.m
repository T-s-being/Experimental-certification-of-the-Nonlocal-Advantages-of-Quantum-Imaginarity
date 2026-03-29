% 2024PRA-Nonlocal Advantages of Quantum Imaginarity

clear all;clc;close all
N_l1_bound = sqrt(5); N_r_bound = 2.02685;
H = [1;0]; V = [0;1];
HH = kron(H,H); VV = kron(V,V); HV = kron(H,V); VH = kron(V,H); 

fai_zheng = (HH+VV)/sqrt(2);
II = kron(eye(2),eye(2));

n = 101;
p = linspace(0,1,n);
best_xall = zeros(n,8);
N_l1 = zeros(1,n);

parfor pi_val = 1:length(p) 
    p_i = p(pi_val);
    rho_AB = p_i*(fai_zheng*fai_zheng') + (1-p_i)/4*II;
    [N_l1(pi_val),best_xall(pi_val,:)] = compute_NAQI(rho_AB);
    fprintf('Progress: %.0f/%.0f\n', pi_val,100); 
end

Experiment_p = [0,0.2,0.4,0.6,0.8,1];
id_p = arrayfun(@(x) find(p==x), Experiment_p);
Experiment_OptimalResult = N_l1(id_p);
Experiment_x = best_xall(id_p,:);
plot(p, N_l1 - N_l1_bound, 'r','HandleVisibility','off'); hold on
scatter(Experiment_p,Experiment_OptimalResult-N_l1_bound,'DisplayName','Experiment dot','MarkerFaceColor','auto')
xlabel('Mixture parameter p');
ylabel('Violation parameter');
legend('show','Location','best')
grid on;
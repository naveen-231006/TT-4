%% =========================================================================
%  TENSEGRITY TUMBLEWEED ROBOT - COMPLETE RESEARCH SIMULATION (FIXED)
%  Design and Analysis of a Tensegrity-Based Tumbleweed Robot
%  with Controlled Locomotion
% =========================================================================
clc; clear; close all;
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');
set(0, 'DefaultAxesZColor', 'k');
set(0, 'DefaultTextColor', 'k');
fprintf('============================================================\n');
fprintf(' Tensegrity Tumbleweed Robot - Research Simulation\n');
fprintf('============================================================\n\n');

%% === SECTION 1: TT-6 GEOMETRY ===
phi_g = (1+sqrt(5))/2;
nodes_raw = [-1,phi_g,0; 1,phi_g,0; -1,-phi_g,0; 1,-phi_g,0;
    0,-1,phi_g; 0,1,phi_g; 0,-1,-phi_g; 0,1,-phi_g;
    phi_g,0,-1; phi_g,0,1; -phi_g,0,-1; -phi_g,0,1];
edge_len = norm(nodes_raw(1,:)-nodes_raw(2,:));
nodes = nodes_raw / edge_len;
N = size(nodes,1);
r_robot = norm(nodes(1,:));
nodes(:,3) = nodes(:,3) + r_robot + 0.05;

bars = [1 4; 2 3; 5 8; 6 7; 9 12; 10 11];
Nb = size(bars,1);

all_e = [];
for i=1:N, for j=i+1:N
    d=norm(nodes(i,:)-nodes(j,:));
    if ~ismember([i j],bars,'rows') && ~ismember([j i],bars,'rows')
        all_e=[all_e; i j d]; end
end, end
thresh = 1.05*min(all_e(:,3));
cables = all_e(all_e(:,3)<=thresh, 1:2);
Nc = size(cables,1);
fprintf('Structure: %d nodes | %d bars | %d cables\n', N, Nb, Nc);

%% === SECTION 2: MAXWELL + FORM-FINDING ===
dof_kin = 3*N-6; members = Nb+Nc; surplus = members-dof_kin;
fprintf('Maxwell: members=%d, DOF=%d, surplus=%d\n', members, dof_kin, surplus);

Cm = zeros(members,N);
for i=1:Nb, Cm(i,bars(i,1))=1; Cm(i,bars(i,2))=-1; end
for i=1:Nc, Cm(Nb+i,cables(i,1))=1; Cm(Nb+i,cables(i,2))=-1; end
q_fd = ones(members,1); q_fd(1:Nb)=-5.0; q_fd(Nb+1:end)=1.5;
D = Cm'*diag(q_fd)*Cm;
res = [norm(D*nodes(:,1)), norm(D*nodes(:,2)), norm(D*nodes(:,3))];
fprintf('Equilibrium residual: %.4f | %.4f | %.4f\n', res);
[~,Ss,~] = svd(D); nz=sum(diag(Ss)<1e-6);
fprintf('Near-zero singular values: %d\n\n', nz);

L0_cables = zeros(Nc,1); L0_bars = zeros(Nb,1);
for i=1:Nc, L0_cables(i)=norm(nodes(cables(i,1),:)-nodes(cables(i,2),:)); end
for i=1:Nb, L0_bars(i)=norm(nodes(bars(i,1),:)-nodes(bars(i,2),:)); end

%% === SECTION 3: MATERIAL PROPERTIES ===
% Cable: Dyneema SK75 (UHMWPE)
d_cable = 0.002; A_cable = pi*(d_cable/2)^2;  % 2mm diameter
E_cable = 100e9; sigma_allow_cable = 900e6; rho_cable = 970;
% Bar: Aluminum 6061-T6 tube (10mm OD, 1mm wall)
d_bar_o = 0.010; d_bar_i = 0.008;
A_bar = pi/4*(d_bar_o^2-d_bar_i^2);
E_bar = 69e9; sigma_yield_bar = 276e6; rho_bar = 2700;

% Mass breakdown
m_bar_each = rho_bar * A_bar * mean(L0_bars);
m_cable_each = rho_cable * A_cable * mean(L0_cables);
m_struct = 6*m_bar_each + 30*m_cable_each;
m_joints = 12 * 0.020;  % 20g per node (3D-printed)
m_total_calc = m_struct + m_joints;
fprintf('--- Materials ---\n');
fprintf('Cable: Dyneema SK75, d=%.0fmm, A=%.4fmm^2, E=%.0fGPa\n', d_cable*1e3, A_cable*1e6, E_cable/1e9);
fprintf('Bar:   Al 6061-T6,  OD=%.0fmm, t=1mm, A=%.4fmm^2, E=%.0fGPa\n', d_bar_o*1e3, A_bar*1e6, E_bar/1e9);
fprintf('Mass breakdown: bars=%.2fkg, cables=%.3fkg, joints=%.2fkg\n', 6*m_bar_each, 30*m_cable_each, m_joints);
fprintf('Total structural mass: %.2f kg (sim: %.2f kg)\n', m_total_calc, 12*0.12);

% Prestrain (0.5% - realistic for Dyneema)
prestrain = 0.005;
L0_cables_pre = L0_cables * (1 - prestrain);
sigma_pre = E_cable * prestrain / 1e6; % MPa
fprintf('Prestress: %.1f MPa (prestrain=%.2f%%)\n\n', sigma_pre, prestrain*100);

%% === SECTION 4: ACTUATION GROUPS ===
cable_groups = zeros(Nc,1);
for i=1:Nc
    mid = 0.5*(nodes(cables(i,1),:)+nodes(cables(i,2),:));
    [~,ax] = max(abs(mid)); cable_groups(i) = ax;
end
gc = histcounts(cable_groups, 0.5:3.5);
fprintf('Actuation groups: X=%d | Y=%d | Z=%d cables\n\n', gc);

%% === SECTION 5: SIMULATION PARAMETERS ===
k_cable = 2000; k_bar = 20000; m_node = 0.12;
dt = 0.05; T_sim = 600; % 30 seconds
rod_radius = 0.03; mu_c = 0.3; mu_v = 0.5;
slope = 0.03; c_damp = 1.5;

% Motor model
wn_mot = 25; zeta_mot = 0.7; u_max = 40;
L_min = 0.50; L_max = 1.40;

% Control rest lengths
L0_ctrl = L0_cables_pre;

%% === SECTION 6: INITIAL STATE ===
% State: [pos(3N); vel(3N); u_motor(Nc); u_dot(Nc)]
x0 = [nodes(:); zeros(3*N,1); zeros(Nc,1); zeros(Nc,1)];

%% === SECTION 7: UKF SETUP ===
nx_u = 9; ny_u = 6;
x_ukf = zeros(nx_u,1); x_ukf(1:3) = mean(nodes,1)';
P_ukf = blkdiag(0.01*eye(3), 0.05*eye(3), 0.02*eye(3));
P_min = 1e-6*eye(nx_u);
Q_ukf = blkdiag(1e-4*eye(3), 5e-4*eye(3), 1e-3*eye(3));
R_ukf = blkdiag(0.004*eye(3), 0.008*eye(2), 0.015);
al=1e-3; be=2; ka=0; lam=al^2*(nx_u+ka)-nx_u;
Wm=[lam/(nx_u+lam), repmat(1/(2*(nx_u+lam)),1,2*nx_u)];
Wc=Wm; Wc(1)=Wc(1)+(1-al^2+be);

%% === SECTION 8: NMPC ===
% Directional steering: target is a far downwind waypoint.
% NMPC steers toward this goal while wind provides propulsion.
Np=6; nu_mpc=3; u_lim=0.15; du_lim=0.06;
Q_mpc = diag([15 50 8 2 2 1 0.1 0.1 0.1]);  % high y-weight for lateral control
R_mpc = 0.3*eye(nu_mpc); Rd_mpc = 1.5*eye(nu_mpc);
x_ref = [50.0; 0; r_robot; 1.4; 0; 0; 0; 0; 0];  % far waypoint + desired velocity
u_prev = zeros(nu_mpc,1);
nmpc_interval = 5; % run NMPC every 5 steps

%% === SECTION 9: SIMULATION LOOP ===
log_com = zeros(T_sim,3); log_vel = zeros(T_sim,3);
log_KE = zeros(T_sim,1); log_PE = zeros(T_sim,1); log_EE = zeros(T_sim,1);
log_u = zeros(T_sim,nu_mpc);
log_bar_f = zeros(T_sim,Nb); log_cable_f = zeros(T_sim,Nc);
log_cov = zeros(T_sim,1); log_shape = zeros(T_sim,3);
log_pow = zeros(T_sim,3); log_Eact = zeros(T_sim,1);
E_act_cum = 0;

x_full = x0;
fprintf('--- Simulation: %d steps x dt=%.3fs = %.1fs ---\n', T_sim, dt, T_sim*dt);

for t = 1:T_sim
    % --- Physics step (ode15s for stiff system) ---
    oopt = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',dt);
    [~,xs] = ode15s(@(tt,xx) dynamics_fn(tt,xx,N,bars,cables,Nb,Nc,...
        k_cable,k_bar,m_node,L0_ctrl,L0_bars,rod_radius,...
        mu_c,mu_v,slope,c_damp,wn_mot,zeta_mot,u_max,L_min,L_max), [0 dt], x_full, oopt);
    x_full = xs(end,:)';
    % Clamp nodes above ground
    pp = reshape(x_full(1:3*N),[N,3]);
    pp(:,3) = max(pp(:,3), 0);
    x_full(1:3*N) = pp(:);

    pos = reshape(x_full(1:3*N),[N,3]);
    vel = reshape(x_full(3*N+1:6*N),[N,3]);
    u_m = x_full(6*N+1:6*N+Nc);

    com = mean(pos,1)'; vel_com = mean(vel,1)';

    % --- Measure cable lengths and forces ---
    clen = zeros(Nc,1); cf = zeros(Nc,1);
    for i=1:Nc
        clen(i) = norm(pos(cables(i,1),:)-pos(cables(i,2),:));
        sl = clen(i)-L0_ctrl(i);
        act = 0.5*(sl+sqrt(sl^2+1e-6));
        cf(i) = max(0, k_cable*act + u_m(i));
    end
    bf = zeros(Nb,1);
    for i=1:Nb
        Lb = norm(pos(bars(i,1),:)-pos(bars(i,2),:));
        bf(i) = k_bar*(Lb-L0_bars(i));
    end

    % Shape deformation per group
    sd = zeros(3,1);
    for g=1:3
        idx = cable_groups==g;
        sd(g) = mean(clen(idx)) - mean(L0_cables(idx));
    end

    % Energy
    KE = 0.5*m_node*sum(vel(:).^2);
    PE = m_node*9.81*sum(pos(:,3));
    EE = 0;
    for i=1:Nc
        EE = EE + 0.5*k_cable*max(0,clen(i)-L0_ctrl(i))^2;
    end

    % Actuator power
    pw = zeros(3,1);
    for g=1:3
        idx=find(cable_groups==g);
        for ii=1:length(idx)
            ci=idx(ii); n1=cables(ci,1); n2=cables(ci,2);
            dv=pos(n2,:)-pos(n1,:); Lc=norm(dv);
            if Lc>1e-9
                vr = dot(vel(n2,:)-vel(n1,:), dv/Lc);
                pw(g) = pw(g) + abs(u_m(ci)*vr);
            end
        end
    end
    E_act_cum = E_act_cum + sum(pw)*dt;

    % --- UKF ---
    if t>1, ac=(com(1)-log_com(t-1,1))/dt-vel_com(1); else, ac=0; end
    z_m = [com+0.003*randn(3,1); vel_com(1:2)+0.008*randn(2,1); ac+0.015*randn];
    P_ukf = max(P_ukf,P_min); P_ukf = 0.5*(P_ukf+P_ukf');
    try Sc=chol((nx_u+lam)*P_ukf,'lower');
    catch, P_ukf=P_ukf+1e-5*eye(nx_u); Sc=chol((nx_u+lam)*P_ukf,'lower'); end
    sp = [x_ukf, bsxfun(@plus,x_ukf,Sc), bsxfun(@minus,x_ukf,Sc)];
    sp_p = zeros(nx_u,2*nx_u+1);
    for k2=1:2*nx_u+1, sp_p(:,k2)=ukf_pred(sp(:,k2),dt); end
    xp = sp_p*Wm'; Pp = Q_ukf;
    for k2=1:2*nx_u+1, d2=sp_p(:,k2)-xp; Pp=Pp+Wc(k2)*(d2*d2'); end
    H = zeros(ny_u,nx_u); H(1:3,1:3)=eye(3); H(4:5,4:5)=eye(2); H(6,4)=1/dt;
    Szz = H*Pp*H'+R_ukf; K_u = Pp*H'/Szz;
    x_ukf = xp + K_u*(z_m-H*xp);
    P_ukf = (eye(nx_u)-K_u*H)*Pp; P_ukf = max(P_ukf,P_min);
    x_ukf(7:9) = sd;

    % --- NMPC (every nmpc_interval steps) ---
    if mod(t,nmpc_interval)==1
        u_opt = nmpc_solve(x_ukf,x_ref,u_prev,Np,nu_mpc,u_lim,du_lim,Q_mpc,R_mpc,Rd_mpc);
        u_prev = u_opt;
    else
        u_opt = u_prev;
    end
    dsc = 0.20;
    for i=1:Nc
        g = cable_groups(i);
        L0_ctrl(i) = L0_cables_pre(i) + dsc*u_opt(g);
        L0_ctrl(i) = max(0.7*L0_cables(i), min(1.15*L0_cables(i), L0_ctrl(i)));
    end

    % --- Log ---
    log_com(t,:)=com'; log_vel(t,:)=vel_com';
    log_KE(t)=KE; log_PE(t)=PE; log_EE(t)=EE;
    log_u(t,:)=u_opt'; log_bar_f(t,:)=bf';
    log_cable_f(t,:)=cf'; log_cov(t)=trace(P_ukf);
    log_shape(t,:)=sd'; log_pow(t,:)=pw'; log_Eact(t)=E_act_cum;

    if mod(t,60)==0
        fprintf('  t=%4d (%.1fs) | CoM=(%.3f,%.3f,%.3f) | KE=%.2f PE=%.2f Eact=%.3f\n',...
            t,t*dt,com(1),com(2),com(3),KE,PE,E_act_cum);
    end
end
fprintf('\nSimulation complete.\n\n');
tv = (1:T_sim)*dt;
pos_fin = reshape(x_full(1:3*N),[N,3]);

%% === SECTION 10: 3D STRUCTURE ===
figure('Name','TT-6 Tensegrity Structure','Position',[50 50 700 600],'Color','w');
hold on; grid on; axis equal;
for i=1:Nb
    p1=pos_fin(bars(i,1),:); p2=pos_fin(bars(i,2),:);
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'r-','LineWidth',4);
end
for i=1:Nc
    p1=pos_fin(cables(i,1),:); p2=pos_fin(cables(i,2),:);
    Lc=norm(p1-p2); strain_c=max(0,(Lc-L0_cables(i))/L0_cables(i));
    sc = min(1, strain_c/0.01);
    col = [sc, 0.3*(1-sc), 1-sc];
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'-','Color',col,'LineWidth',1.5);
end
scatter3(pos_fin(:,1),pos_fin(:,2),pos_fin(:,3),80,'k','filled');
% Ground plane
xl=xlim; yl=ylim;
patch([xl(1) xl(2) xl(2) xl(1)],[yl(1) yl(1) yl(2) yl(2)],[0 0 0 0],...
    [0.8 0.9 0.8],'FaceAlpha',0.3,'EdgeColor','none');
view(45,30);
title('TT-6 Tensegrity - Final Configuration');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
legend('Bars (Al 6061-T6)','','','','','','Cables (Dyneema SK75)','Location','best');
set(gcf, 'Color', 'w');
saveas(gcf,'fig_structure_3d.png');

%% === SECTION 11: SYSTEM ANALYSIS (6-panel) ===
figure('Name','System Analysis','Position',[30 30 1400 900],'Color','w');

subplot(2,3,1);
plot(log_com(:,1),log_com(:,2),'b-','LineWidth',2); hold on;
plot(log_com(1,1),log_com(1,2),'go','MarkerSize',10,'LineWidth',2);
plot(x_ref(1),x_ref(2),'rp','MarkerSize',14,'LineWidth',2);
xlabel('X [m]'); ylabel('Y [m]');
title('Center of Mass Trajectory (XY Plane)');
legend('Path','Start','Target','Location','best'); grid on; axis equal;

subplot(2,3,2);
plot(tv,log_com(:,1),'b',tv,log_com(:,2),'r',tv,log_com(:,3),'g','LineWidth',1.5);
yline(x_ref(1),'b--','x_{ref}');
xlabel('Time [s]'); ylabel('Position [m]');
title('CoM Position Components'); legend('x','y','z'); grid on;

subplot(2,3,3);
err_xy = vecnorm(log_com(:,1:2)-x_ref(1:2)',2,2);
plot(tv, err_xy, 'm-', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Error [m]');
title('Tracking Error (XY distance to target)'); grid on;

subplot(2,3,4);
E_mech = log_KE + log_PE + log_EE;
plot(tv,log_KE,'b',tv,log_PE,'r',tv,log_EE,'g',tv,E_mech,'k--','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Energy [J]');
title('Energy Components'); legend('KE','PE_{grav}','PE_{elastic}','E_{mech}'); grid on;

subplot(2,3,5);
plot(tv,log_u(:,1),'b',tv,log_u(:,2),'r',tv,log_u(:,3),'g','LineWidth',1.5);
yline(u_lim,'k--'); yline(-u_lim,'k--');
xlabel('Time [s]'); ylabel('u [m]');
title('NMPC Control Inputs'); legend('X-group','Y-group','Z-group'); grid on;

subplot(2,3,6);
semilogy(tv,log_cov,'m','LineWidth',1.5);
xlabel('Time [s]'); ylabel('tr(P)');
title('UKF Covariance Trace (State Estimation)'); grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_system_analysis.png');

%% === SECTION 12: STRESS & UTILISATION ===
% Physical stress from simulation forces
% Cable: prestress + dynamic (quasi-static force / area)
sigma_cable_MPa = mean(log_cable_f,1)' / A_cable / 1e6 + sigma_pre;
% Bar: includes compression from cable pretension
% Each bar end connects to cables pulling inward -> net axial compression
F_pre_cable = sigma_pre * 1e6 * A_cable;  % prestress force per cable [N]
sigma_bar_comp = F_pre_cable * 4 / A_bar / 1e6;  % ~4 cables load each bar
sigma_bar_MPa = abs(mean(log_bar_f,1))' / A_bar / 1e6 + sigma_bar_comp;
util_cable = sigma_cable_MPa / (sigma_allow_cable/1e6) * 100;
util_bar = sigma_bar_MPa / (sigma_yield_bar/1e6) * 100;

figure('Name','Structural Stress Analysis','Position',[60 60 1200 500],'Color','w');
subplot(1,3,1);
bar(sigma_cable_MPa,'FaceColor',[0.2 0.5 0.9]);
yline(sigma_allow_cable/1e6,'r--','Allowable','LineWidth',1.5);
xlabel('Cable Index'); ylabel('Stress [MPa]');
title(sprintf('Cable Stress (Dyneema SK75)\nPrestress = %.0f MPa',sigma_pre)); grid on;

subplot(1,3,2);
barh(sigma_bar_MPa,'FaceColor',[0.9 0.3 0.2]);
xline(sigma_yield_bar/1e6,'r--','Yield','LineWidth',1.5);
xlabel('Stress [MPa]'); ylabel('Bar Index');
title('Bar Stress (Al 6061-T6)'); grid on;

subplot(1,3,3);
bar([max(util_cable); max(util_bar)],'FaceColor',[0.3 0.7 0.4]);
xline(100,'r--','100% Yield','LineWidth',1.5);
set(gca,'XTickLabel',{'Cable max','Bar max'});
ylabel('Utilisation [%]');
title(sprintf('Structural Utilisation\nCable: %.1f%% | Bar: %.1f%%',max(util_cable),max(util_bar)));
grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_stress_analysis.png');

%% === SECTION 13: ENERGY & POWER ===
figure('Name','Energy and Actuation','Position',[50 50 1100 450],'Color','w');
subplot(1,2,1);
plot(tv,E_mech,'b',tv,log_Eact,'r','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Energy [J]');
title('Mechanical Energy vs Actuator Input');
legend('E_{mech}','E_{actuator} (cumulative)'); grid on;

subplot(1,2,2);
plot(tv,log_pow(:,1),'b',tv,log_pow(:,2),'r',tv,log_pow(:,3),'g','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Power [W]');
title('Instantaneous Actuator Power'); legend('X-group','Y-group','Z-group'); grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_energy_power.png');

%% === SECTION 14: SHAPE DEFORMATION ===
figure('Name','Shape Mode Deformation','Position',[70 70 800 400],'Color','w');
plot(tv,log_shape(:,1)*1e3,'b',tv,log_shape(:,2)*1e3,'r',tv,log_shape(:,3)*1e3,'g','LineWidth',1.5);
xlabel('Time [s]'); ylabel('\Delta L [mm]');
title('Cable Group Mean Deformation'); legend('X-group','Y-group','Z-group'); grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_shape_deformation.png');

%% === SECTION 15: MATERIAL COMPARISON ===
mat_names = {'Dyneema SK75','Kevlar 49','Steel Wire','Carbon Fiber'};
mat_E = [100, 112, 200, 230]; % GPa
mat_sigma = [3400, 3600, 1800, 3500]; % MPa UTS
mat_rho = [970, 1440, 7800, 1750]; % kg/m3
spec_str = mat_sigma ./ mat_rho; % kN*m/kg

figure('Name','Material Comparison','Position',[80 80 900 500],'Color','w');
subplot(1,2,1);
bar(mat_E,'FaceColor',[0.3 0.5 0.8]);
set(gca,'XTickLabel',mat_names);
ylabel("Young's Modulus [GPa]"); title('Cable Material Stiffness'); grid on;

subplot(1,2,2);
bar(spec_str,'FaceColor',[0.8 0.4 0.2]);
set(gca,'XTickLabel',mat_names);
ylabel('Specific Strength [kN\cdotm/kg]');
title('Cable Material Specific Strength'); grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_material_comparison.png');

%% === SECTION 16: PHASE PORTRAIT ===
figure('Name','Phase Portrait','Position',[90 90 600 500],'Color','w');
plot(log_com(:,1),log_vel(:,1),'b','LineWidth',1.5); hold on;
plot(log_com(:,2),log_vel(:,2),'r','LineWidth',1.5);
xlabel('Position [m]'); ylabel('Velocity [m/s]');
title('Phase Portrait (CoM)'); legend('X','Y'); grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_phase_portrait.png');

%% === SECTION 17: PERFORMANCE SUMMARY ===
final_err = norm(log_com(end,1:2)'-x_ref(1:2));
dist_trav = sum(vecnorm(diff(log_com(:,1:2)),2,2));
ci = find(err_xy<0.1,1); ct = NaN; if ~isempty(ci), ct=ci*dt; end

fprintf('====================== PERFORMANCE SUMMARY ======================\n');
fprintf('  Simulation time:       %.1f s (%d steps)\n', T_sim*dt, T_sim);
fprintf('  Final CoM (x,y,z):     %.4f  %.4f  %.4f m\n', log_com(end,:));
fprintf('  Target (x,y):          %.4f  %.4f m\n', x_ref(1), x_ref(2));
fprintf('  Final tracking error:  %.4f m\n', final_err);
fprintf('  Convergence time:      '); if isnan(ct), fprintf('Not converged\n'); else, fprintf('%.2f s\n',ct); end
fprintf('  Distance travelled:    %.4f m\n', dist_trav);
fprintf('  Mean speed:            %.4f m/s\n', mean(vecnorm(log_vel(:,1:2),2,2)));
fprintf('  Total actuator energy: %.4f J\n', E_act_cum);
fprintf('  Locomotion efficiency: %.4f m/J\n', dist_trav/max(E_act_cum,1e-9));
fprintf('  Cable prestress:       %.1f MPa\n', sigma_pre);
fprintf('  Max cable stress:      %.1f MPa (util: %.1f%%)\n', max(sigma_cable_MPa), max(util_cable));
fprintf('  Max bar stress:        %.1f MPa (util: %.1f%%)\n', max(sigma_bar_MPa), max(util_bar));
fprintf('  Cable safety factor:   %.2f\n', (sigma_allow_cable/1e6)/max(sigma_cable_MPa));
fprintf('  Bar safety factor:     %.2f\n', (sigma_yield_bar/1e6)/max(sigma_bar_MPa));
fprintf('  Energy drift:          %.4f J (%.1f%%)\n', abs(E_mech(end)-E_mech(1)), abs(E_mech(end)-E_mech(1))/E_mech(1)*100);
fprintf('  Material (cables):     Dyneema SK75 (UHMWPE)\n');
fprintf('  Material (bars):       Aluminum 6061-T6\n');
fprintf('================================================================\n');
fprintf('Figures saved.\n');

%% === SECTION 18: BASELINE COMPARISON (NO CONTROLLER) ===
fprintf('\n--- Running baseline (no NMPC) for comparison ---\n');
x_base = x0; L0_base = L0_cables_pre;
log_com_base = zeros(T_sim,3);
for t = 1:T_sim
    oopt = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',dt);
    [~,xs] = ode15s(@(tt,xx) dynamics_fn(tt,xx,N,bars,cables,Nb,Nc,...
        k_cable,k_bar,m_node,L0_base,L0_bars,rod_radius,...
        mu_c,mu_v,slope,c_damp,wn_mot,zeta_mot,u_max,L_min,L_max), [0 dt], x_base, oopt);
    x_base = xs(end,:)';
    pp = reshape(x_base(1:3*N),[N,3]);
    pp(:,3) = max(pp(:,3), 0);
    x_base(1:3*N) = pp(:);
    log_com_base(t,:) = mean(pp,1);
    if mod(t,120)==0
        fprintf('  Baseline t=%4d | CoM=(%.3f,%.3f,%.3f)\n',t,log_com_base(t,:));
    end
end
base_lat = log_com_base(end,2);
ctrl_lat = log_com(end,2);
fprintf('  Baseline lateral drift:    %.3f m\n', base_lat);
fprintf('  Controlled lateral drift:  %.3f m\n', ctrl_lat);
fprintf('  Lateral reduction:         %.1f%%\n', (1-abs(ctrl_lat)/max(abs(base_lat),1e-6))*100);

figure('Name','Baseline Comparison','Position',[50 50 1000 450],'Color','w');
subplot(1,2,1);
plot(log_com_base(:,1),log_com_base(:,2),'r--','LineWidth',1.5); hold on;
plot(log_com(:,1),log_com(:,2),'b-','LineWidth',2);
plot(log_com(1,1),log_com(1,2),'go','MarkerSize',10,'LineWidth',2);
xlabel('X [m]'); ylabel('Y [m]');
title('Trajectory: Controlled vs Passive');
legend('Passive (no NMPC)','Controlled (NMPC)','Start','Location','best');
grid on; axis equal;
subplot(1,2,2);
lat_base = log_com_base(:,2);
lat_ctrl = log_com(:,2);
plot(tv, lat_base, 'r--', tv, lat_ctrl, 'b-', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('Lateral Drift Y [m]');
title('Lateral Deviation Over Time');
legend('Passive','Controlled'); grid on;
set(gcf, 'Color', 'w');
saveas(gcf,'fig_baseline_comparison.png');
fprintf('Baseline comparison figure saved.\n');

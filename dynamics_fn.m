function dxdt = dynamics_fn(~, x, N, bars, cables, Nb, Nc, ...
    k_cable, k_bar, m_node, L0_ctrl, L0_bars, rod_radius, ...
    mu_c, mu_v, slope, c_damp, wn_mot, zeta_mot, u_max, L_min, L_max)
% Tensegrity dynamics: cables + bars + contact + friction + motor

pos = reshape(x(1:3*N), [N,3]);
vel = reshape(x(3*N+1:6*N), [N,3]);
u_m = x(6*N+1:6*N+Nc);
u_d = x(6*N+Nc+1:end);

F = zeros(N,3);

% Cable forces (tension-only with softplus activation)
for i = 1:Nc
    n1=cables(i,1); n2=cables(i,2);
    dvec = pos(n2,:)-pos(n1,:);
    L = norm(dvec);
    if L < 1e-9, continue; end
    dir = dvec/L;
    slack = L - L0_ctrl(i);
    act = 0.5*(slack + sqrt(slack^2 + 1e-6));
    f_mag = max(0, k_cable*act + u_m(i));
    if L < L_min, f_mag = 0; end
    if L > L_max, f_mag = u_max; end
    f = f_mag * dir;
    F(n1,:) = F(n1,:) + f;
    F(n2,:) = F(n2,:) - f;
end

% Bar rigidity (bilateral penalty - high stiffness)
for i = 1:Nb
    n1=bars(i,1); n2=bars(i,2);
    dvec = pos(n2,:)-pos(n1,:);
    L = norm(dvec);
    if L < 1e-9, continue; end
    f = (k_bar*(L-L0_bars(i))/L) * dvec;
    F(n1,:) = F(n1,:) + f;
    F(n2,:) = F(n2,:) - f;
    % Bar damping
    vrel = vel(n2,:)-vel(n1,:);
    fd = 5.0 * dot(vrel, dvec/L) * (dvec/L);
    F(n1,:) = F(n1,:) + fd;
    F(n2,:) = F(n2,:) - fd;
end

% Rod rotational inertia
for i = 1:Nb
    n1=bars(i,1); n2=bars(i,2);
    rv = pos(n2,:)-pos(n1,:); Lr = norm(rv);
    if Lr < 1e-9, continue; end
    rh = rv/2;
    vr = vel(n2,:)-vel(n1,:);
    omega = cross(rh, vr)/(dot(rh,rh)+1e-9);
    m_rod = m_node*2; I_rod = (m_rod*Lr^2)/12;
    Iw = I_rod*omega;
    tau = -cross(omega, Iw);
    fr = cross(tau, rh)/(dot(rh,rh)+1e-9);
    F(n1,:) = F(n1,:) - 0.05*fr;
    F(n2,:) = F(n2,:) + 0.05*fr;
end

% Torsional stiffness
k_tor = 2.0;
for i = 1:Nb
    for j = i+1:Nb
        a1=bars(i,1); a2=bars(i,2); b1=bars(j,1); b2=bars(j,2);
        sh = any(ismember(cables(:,1),[a1 a2])&ismember(cables(:,2),[b1 b2])) || ...
             any(ismember(cables(:,1),[b1 b2])&ismember(cables(:,2),[a1 a2]));
        if sh
            ui=(pos(a2,:)-pos(a1,:)); ui=ui/(norm(ui)+1e-9);
            uj=(pos(b2,:)-pos(b1,:)); uj=uj/(norm(uj)+1e-9);
            ct=max(-1,min(1,dot(ui,uj)));
            th=acos(ct); cr=cross(ui,uj); nc=norm(cr);
            if nc>1e-6
                ax=cr/nc; tq=k_tor*th;
                F(a1,:)=F(a1,:)-0.01*tq*ax; F(a2,:)=F(a2,:)+0.01*tq*ax;
                F(b1,:)=F(b1,:)+0.01*tq*ax; F(b2,:)=F(b2,:)-0.01*tq*ax;
            end
        end
    end
end

% Rod-rod collision
k_col = 800;
for i = 1:Nb
    for j = i+1:Nb
        a1=bars(i,1);a2=bars(i,2);b1=bars(j,1);b2=bars(j,2);
        [d,cd] = seg_dist_fn(pos(a1,:),pos(a2,:),pos(b1,:),pos(b2,:));
        if d < 2*rod_radius && d > 1e-6
            pen = 2*rod_radius-d;
            fc = k_col*pen*cd;
            F(a1,:)=F(a1,:)+fc/2; F(a2,:)=F(a2,:)+fc/2;
            F(b1,:)=F(b1,:)-fc/2; F(b2,:)=F(b2,:)-fc/2;
        end
    end
end

% Gravity on slope
g_vec = [9.81*sin(slope), 0, -9.81*cos(slope)];
F = F + m_node*repmat(g_vec, N, 1);

% Ground contact with Coulomb friction
k_gnd = 5000; d_gnd = 150;
for i = 1:N
    gz = slope*pos(i,1);
    if pos(i,3) < gz
        pen = gz - pos(i,3);
        Fn = max(0, k_gnd*pen - d_gnd*vel(i,3));
        F(i,3) = F(i,3) + Fn;
        vt = vel(i,1:2); vtm = norm(vt)+1e-6;
        Ft = -mu_c*Fn*(vt/vtm) - mu_v*vt;
        F(i,1:2) = F(i,1:2) + Ft;
    end
end

% Rolling resistance
cv = mean(vel,1); vr = norm(cv(1:2));
if vr > 1e-4
    W = N*m_node*9.81;
    fr = 0.005*W/N;  % reduced rolling resistance
    for i = 1:N
        F(i,1) = F(i,1) - fr*(cv(1)/vr);
        F(i,2) = F(i,2) - fr*(cv(2)/vr);
    end
end

% Aerodynamic wind force (tumbleweed driving force)
% Wind applies drag to exposed (upper) nodes
v_wind = [8.0, 1.0, 0]; % wind velocity [m/s]
Cd = 1.0; rho_air = 1.225; A_node = 0.07; % drag area per node (~0.84m2 total)
com_z = mean(pos(:,3));
for i = 1:N
    if pos(i,3) > com_z - 0.1  % exposed to wind
        v_rel = v_wind - vel(i,:);
        F_drag = 0.5*rho_air*Cd*A_node*norm(v_rel)*v_rel;
        F(i,:) = F(i,:) + F_drag;
    end
end

% Structural damping (light)
F = F - 0.5*vel;

% Motor dynamics (2nd order, PD control)
Kp = 80; Kd = 12;
du_m = zeros(Nc,1);
du_d = zeros(Nc,1);
for i = 1:Nc
    n1=cables(i,1); n2=cables(i,2);
    Lc = norm(pos(n2,:)-pos(n1,:));
    ep = L0_ctrl(i) - Lc;
    dep = -dot(vel(n2,:)-vel(n1,:), (pos(n2,:)-pos(n1,:))/(Lc+1e-9));
    uc = Kp*ep + Kd*dep;
    uc = max(min(uc, u_max), -u_max);
    du_m(i) = u_d(i);
    du_d(i) = wn_mot^2*(uc-u_m(i)) - 2*zeta_mot*wn_mot*u_d(i);
end

acc = F/m_node;
dxdt = [vel(:); acc(:); du_m; du_d];
end

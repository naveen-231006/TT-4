function u_opt = nmpc_solve(x0, x_ref, u_prev, Np, nu, u_lim, du_lim, Q, R, Rd)
% NMPC solver using SQP (fmincon)
u0_v = repmat(u_prev, Np, 1);
lb_v = -u_lim*ones(nu*Np,1);
ub_v =  u_lim*ones(nu*Np,1);

% Rate constraints
In = eye(nu);
Ad = zeros(nu*(Np-1), nu*Np);
for k = 1:Np-1
    Ad((k-1)*nu+1:k*nu, (k-1)*nu+1:k*nu) = -In;
    Ad((k-1)*nu+1:k*nu, k*nu+1:(k+1)*nu) =  In;
end
bd = du_lim*ones(nu*(Np-1),1);

opts = optimoptions('fmincon','Display','off','MaxIterations',25,...
    'Algorithm','sqp','OptimalityTolerance',1e-4);
u_seq = fmincon(@(u)nmpc_cost(u,x0,x_ref,u_prev,Np,nu,Q,R,Rd), u0_v,...
    [Ad;-Ad],[bd;bd],[],[],lb_v,ub_v,[],opts);
u_opt = u_seq(1:nu);
end

function J = nmpc_cost(u_seq, x0, x_ref, u_prev_in, Np, nu, Q, R, Rd)
x = x0; J = 0; upk = u_prev_in;
dt_p = 0.05;
k_eff = 80.0; mass_tot = 12*0.12; drag = 0.25;
% Wind feedforward: approximate net wind acceleration on exposed nodes
% v_wind=[8,1,0], ~8 of 12 nodes exposed, Cd=1, A_node=0.07, rho=1.225
a_wind = [0.28; 0.04; 0];  % estimated wind acceleration [m/s^2]
for k = 1:Np
    uk = u_seq((k-1)*nu+1:k*nu);
    % Predict next state (reduced-order with wind feedforward)
    a_ctrl = (k_eff/mass_tot)*uk;
    x_n = x;
    x_n(1:3) = x(1:3) + dt_p*x(4:6);
    x_n(4) = x(4) + dt_p*(a_ctrl(1) + a_wind(1) - drag*x(4));
    x_n(5) = x(5) + dt_p*(a_ctrl(2) + a_wind(2) - drag*x(5));
    x_n(6) = max(-0.5, x(6)+dt_p*(a_ctrl(3)-9.81-drag*x(6)));
    tau_s = 5.0;
    x_n(7) = x(7)*(1-dt_p/tau_s)+0.01*uk(1);
    x_n(8) = x(8)*(1-dt_p/tau_s)+0.01*uk(2);
    x_n(9) = x(9)*(1-dt_p/tau_s)+0.01*uk(3);
    x = x_n;
    % Terminal weight
    w = 1 + double(k==Np)*2;
    e = x - x_ref;
    J = J + w*(e'*Q*e) + uk'*R*uk + (uk-upk)'*Rd*(uk-upk);
    upk = uk;
end
end

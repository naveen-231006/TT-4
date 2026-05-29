function x_next = ukf_pred(x_in, dt_in)
% UKF prediction model for 9D state
x_next = x_in;
x_next(1:3) = x_in(1:3) + dt_in*x_in(4:6);
drag = 0.25;
x_next(4) = x_in(4) - dt_in*drag*x_in(4);
x_next(5) = x_in(5) - dt_in*drag*x_in(5);
x_next(6) = max(-0.5, x_in(6) - dt_in*(9.81+drag*x_in(6)));
tau_s = 5.0;
x_next(7:9) = x_in(7:9)*(1-dt_in/tau_s);
end

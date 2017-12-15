function heuns(a, x_init)
    t_init = 0; step_size = 0.01; t_final = 0.05;
%     x_, y_ = heun's method, y2_ = euler's method
    [u1,v1] = solve_system(a, x_init);
    
    [t_, x1_, x2_, xe1_, xe2_] = heuns_method(x_init, t_init, step_size, t_final, u1, v1);
    t_actual = (t_init : 0.0001 : t_final);
    xsol1 = -0.02 * exp(-25 * t_actual) + 0.04;
    xsol2 = 0.02 * exp(-25 * t_actual) + 0.06;
    plot(t_, xe1_, t_, xe2_, t_actual, xsol1, t_actual, xsol2, t_, x1_, t_, x2_), legend("x1 eulers approx", "x2 eulers approx", "x1 solution", "x2 solution", "x2 heun's approx", "x1 heun's approx");
end

function [t, x1, x2, xe1, xe2] = heuns_method(x_init, t_init, step_size, t_final, u, v)
    n = 1; % counter
    t(n) = t_init;
    xe1(n) = x_init(1);
    xe2(n) = x_init(2);
    xh1(n) = x_init(1);
    xh2(n) = x_init(2);
    x1(n) = x_init(1);
    x2(n) = x_init(2);
    while (t(n) < t_final)
        n = n + 1;
        t(n) = t(n - 1) + step_size;
        xe1(n) = xe1(n - 1) + step_size * diffeq_x1(xe1(n - 1), xe2(n - 1), t(n - 1)); % euler's method
        xe2(n) = xe2(n - 1) + step_size * diffeq_x2(xe1(n - 1), xe2(n - 1), t(n - 1)); % euler's method
        xh1(n) = x1(n - 1) + step_size * diffeq_x1(x1(n - 1), x2(n - 1), t(n - 1));
        xh2(n) = x2(n - 1) + step_size * diffeq_x2(x1(n - 1), x2(n - 1), t(n - 1));
        x1(n) = x1(n - 1) + step_size / 2 * (diffeq_x1(x1(n - 1), x2(n - 1), t(n - 1)) + diffeq_x1(xh1(n), xh2(n), t(n))); % euler's method
        x2(n) = x2(n - 1) + step_size / 2 * (diffeq_x2(x1(n - 1), x2(n - 1), t(n - 1)) + diffeq_x2(xh1(n), xh2(n), t(n))); % euler's method
        error_heuns1 = abs(eval(u(t(n))) - x1(n)) * 100 / eval(u(t(n)));
        error_heuns2 = abs(eval(v(t(n))) - x2(n)) * 100 / eval(v(t(n)));
        error_eulers1 = abs(eval(u(t(n))) - xe1(n)) * 100 / eval(u(t(n)));
        error_eulers2 = abs(eval(v(t(n))) - xe2(n)) * 100 / eval(v(t(n)));
        disp("error heuns x1: " + error_heuns1 + " %; error heuns x2: " + error_heuns2 + " %");
        disp("error eulers x1: " + error_eulers1 + " %; error eulers x2: " + error_eulers2 + " %");
        disp("heun's approximation is more accurate than eulers by " + (error_eulers1 - error_heuns1) + "%");
        disp(" ");
    end
end

function [u, v] = solve_system(a, x_init)
    syms u(t) v(t);

    ode1 = diff(u) == a(1,1)*u + a(1,2)*v;
    ode2 = diff(v) == a(2,1)*u + a(2,2)*v;
    odes = [ode1; ode2];
    
    [uSol(t), vSol(t)] = dsolve(odes);
    
    cond1 = u(0) == x_init(1);
    cond2 = v(0) == x_init(2);
    conds = [cond1; cond2];
    [uSol(t), vSol(t)] = dsolve(odes,conds);
    
    fplot(uSol, [0,3])
    hold on
    fplot(vSol, [0,3])
    grid on
    legend('uSol','vSol','Location','best')
    
    u = uSol;
    v = vSol;
end

% first order differential equation
function x1_prime = diffeq_x1(x1, x2, t)
    x1_prime = -15 * x1 + 10 * x2;
end

function x2_prime = diffeq_x2(x1, x2, t)
    x2_prime = 15 * x1 - 10 * x2;
end
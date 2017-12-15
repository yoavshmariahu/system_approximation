function heuns_op(a, x_init, t_init, t_final, step_size)
% a contains n by n matrix of coefficients of the system of equations
% x_init contains n size vector of initial values for system of equations
% approximates system from t_init to t_final, incrementing by step_size
    t_ = t_init : step_size : t_final;
    [X_] = approx_system(a, x_init, t_, step_size);
    plot(t_, X_);
end

function [X] = approx_system(a, x_init, t_, step_size)
% matrix X contains optomized Heun's approximation values
% matrix Xh contains helper approximation values for X
% matrix Xe contains euler's method approximation values
    X = zeros(1, length(t_));
    Xh = zeros(1, length(t_));
    Xe = zeros(1, length(t_));
    for n = 1:length(a)
        X(n, 1) = x_init(n);
        Xh(n, 1) = x_init(n);
        Xe(n, 1) = x_init(n);
    end
    for n = 2:length(t_)
        for m = 1:length(a)
            Xe(m, n) = Xe(m, n - 1) + step_size * (a(m,:) * Xe(:,n - 1));
            Xh(m, n) = X(m, n - 1) + step_size * (a(m,:) * Xh(:,n - 1));
        end
        for m = 1:length(a)
            X(m, n) = X(m, n - 1) + step_size / 2 * ((a(m,:) * X(:,n - 1)) + (a(m,:) * Xh(:,n)));
        end
    end
end

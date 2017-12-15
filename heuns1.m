function heuns()
    x_init = 0; y_init = 0; step_size = 0.5; x_final = 3;
    % x_, y_ = heun's method, y2_ = euler's method
    [x_, y_, y2_] = heuns_method(x_init, y_init, step_size, x_final);
    a = (x_init : 0.0001 : x_final);
    b = a.^3;
    plot(x_, y_, x_, y2_), legend("heuns approximation", "euler's approximation", "y = x^2");
end

function [x, y, y2] = heuns_method(x_init, y_init, step_size, x_final)
    n = 1; % counter
    x(n) = x_init;
    y(n) = y_init;
    y2(n) = y_init;
    while (x(n) < x_final)
        n = n + 1;
        x(n) = x(n - 1) + step_size;
        y2(n) = y(n - 1) + step_size * diffeq(x(n - 1), y(n - 1)); % euler's method
        y(n) = y(n - 1) + step_size / 2 * (diffeq(x(n - 1), y(n - 1)) + diffeq(x(n), y2(n))); % improved euler's method
        disp("x: " + x(n) + "; euler's: y = " + y2(n) + ", heun's: y = " + y(n) + " actual: y = " + eq(x(n), y_init));
    end
end

% first order differential equation
function a = diffeq(x, y)
    a = 3 * x.^2;
end

% equation
function b = eq(x, y_init)
    b = x.^3 + y_init;
end
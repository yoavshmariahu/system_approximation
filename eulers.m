function eulers()
[x_,y_] = heuns_method(0,0,0.3,1);
a = (0:0.001:4);
b = a.^2;
plot(x_, y_, a, b), legend('heuns approx', 'y = x^2');
end

function [x,y] = heuns_method(x_init, y_init, step_size, n)
    x(n) = x_init;
    y(n) = diffeq(x_init, y_init) * step_size + y_init;
    while (x(n) < 4)
        n = n + 1;
        x(n) = x(n - 1) + step_size;
        r = y(n - 1) + step_size * diffeq(x(n - 1), y(n - 1));
        y(n) = y(n - 1) + step_size / 2 * (diffeq(x(n - 1), y(n - 1)) + diffeq(x(n), r));
        disp("n : " + n + " (" + x(n) + ", " + y(n) + ")");
    end
end

function a = diffeq(x, y)
a = 2 * x;
end
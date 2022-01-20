clear
% Define the function
f = @(x) x.^3 - 3*x +1;
% Plot
x = linspace(-3, 3, 100);
plot(x, f(x), [-3 3], [0, 0])
% Check analytical roots
f(-2*cos(pi/9))
f(-sqrt(3)*sin(pi/9) + cos(pi/9))
% Methods for finding the roots`
%  - Brent's methods
fzero(f, -3)
fzero(f, 0)
fzero(f, 2)
%  - Newton's method
df = @(x) 3*x.^2 - 3;
x(1) = -3.0;
err_tol  = 0.00001;
max_iter = 10;
for i = 1:max_iter
    x(i+1) = x(i)-((f(x(i)))/df(x(i)));
    err(i) = abs((x(i+1)-x(i))/x(i));
    if err(i) < err_tol
        break
    end
end
root=x(i);
disp(root)
%  - Fixed-point iteration
%    Note: we need to multiply by a factor to satisfy the convergence 
%          condition, |g'| < 1, for fixed-point iteration 
g = @(x) 0.1*((x.^3)/3 + 1/3);
x = linspace(-3, 3, 100);
figure
plot(x,g(x), x, x/10)
%
y0 = -3./10;
err_tol  = 0.00001;
max_iter = 100;
for i = 1:max_iter
    y = 10*g(y0);
    err = abs(y-y0);
    if err < err_tol
        break
    end
    y0 = y;
end
disp(y);

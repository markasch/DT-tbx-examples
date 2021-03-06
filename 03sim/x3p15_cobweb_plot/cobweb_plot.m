function cobweb_plot(f,a,b,x0,N)
% Generate the cobweb plot associated with
% the orbits of x_n+1=f(x_n).
% Input:
%    N is the number of iterates
%    (a,b) is the interval
%    x0 is the initial point
%    @f to pass the function 
%
% Generate N linearly spaced values on (a,b)
x = linspace(a,b,N);  
% which we use to plot the function y=f(x)
y = f(x);
% Turn hold on to gather up all plots in one
hold on;
plot(x,y);     % plot the function
plot(x,x,'r'); % plot the diagonal line
x(1) = x0;     % plot orbit starting at x0
for i = 1:N
   x(i+1) = f(x(i));
   line([x(i),x(i)],[x(i),x(i+1)]);
   line([x(i),x(i+1)],[x(i+1),x(i+1)]);
end
hold off;
% Function f.m
function ret=f(x) 
r   = 2; % try other values
ret = r*x.*(1-x);
end
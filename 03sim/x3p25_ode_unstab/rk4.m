function [ t, y ] = rk4 ( yprime, tspan, y0, n )

%*****************************************************************************80
%
%% rk4 uses an explicit Runge-Kutta scheme of order 4.
%
%  Discussion:
%
%    This function approximates the solution of a differential
%    equation of the form:
%      y' = yprime(t,y)
%      y(tspan(1)) = y0
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 March 2020
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    function handle yprime: defines the derivative function of the form
%      value = yprime ( t, y )
%
%    real tspan(2): contains the initial and final times.
%
%    real y0(m): contains the initial condition.
%
%    integer n: the number of steps to take.
%
%  Output:
%
%    real t(n+1): the times.
%
%    real y(n+1,m): the estimated solution values.
%
  m = length ( y0 );
  t = zeros ( n + 1, 1 );
  y = zeros ( n + 1, m );

  dt = ( tspan(2) - tspan(1) ) / n;

  t(1,1) = tspan(1);
  y(1,:) = y0(:);

  a = [       0.0,          0.0,          0.0,        0.0;
             0.25,          0.0,          0.0,        0.0;
         4.0/81.0,    32.0/81.0,          0.0,        0.0;
        57.0/98.0, -432.0/343.0, 1053.0/686.0,        0.0;
          1.0/6.0,          0.0,    27.0/52.0, 49.0/156.0];

  b = [ 43.0/288.0, 0.0, 243.0/416.0, 343.0/1872.0, 1.0/12.0 ];
  c = [ 0.0, 0.25, 4.0/9.0, 6.0/7.0, 1.0 ];
 
  for i = 1 : n

    k1 = dt * yprime ( t(i) + c(1) * dt, y(i,:) );
    k2 = dt * yprime ( t(i) + c(2) * dt, y(i,:) + a(2,1) * k1' );
    k3 = dt * yprime ( t(i) + c(3) * dt, y(i,:) + a(3,1) * k1' ...
      + a(3,2) * k2' );
    k4 = dt * yprime ( t(i) + c(4) * dt, y(i,:) + a(4,1) * k1' ...
      + a(4,2) * k2' + a(4,3) * k3' );
    k5 = dt * yprime ( t(i) + c(5) * dt, y(i,:) + a(5,1) * k1' ...
      + a(5,2) * k2' + a(5,3) * k3' + a(5,4) * k4' );

    y4 = y(i,:) + b(1) * k1' + b(2) * k2' + b(3) * k3' + b(4) * k4' + b(5) * k5';

    t(i+1,1) = t(i,1) + dt;
    y(i+1,:) = y4;

  end

  return
end


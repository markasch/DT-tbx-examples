function [t,u]=beuler(odefun,tspan,y0,Nh,varargin)
%BEULER Solves differential equations using the
%  backward Euler method.
%  [T,Y]=BEULER(ODEFUN,TSPAN,Y0,NH) with TSPAN=[T0,TF]
%  integrates the system of differential equations
%  y'=f(t,y) from time T0 to TF with initial condition
%  Y0 using the backward Euler method on an equispaced
%  grid of NH intervals.
%  Function ODEFUN(T,Y) must return a vector, whose
%  elements hold the evaluation of f(t,y), of the
%  same dimension of Y.
%  Each row in the solution array Y corresponds to a
%  time returned in the  column vector T.
%  [T,Y] = BEULER(ODEFUN,TSPAN,Y0,NH,P1,P2,...) passes
%  the additional parameters P1,P2,... to the function
%  ODEFUN as ODEFUN(T,Y,P1,P2...).
tt=linspace(tspan(1),tspan(2),Nh+1);
y=y0(:); % always create a vector column
u=y.';
global glob_h glob_t glob_y glob_odefun;
glob_h=(tspan(2)-tspan(1))/Nh;
glob_y=y;
glob_odefun=odefun;
glob_t=tt(2);

if ( exist('OCTAVE_VERSION') )
o_ver=OCTAVE_VERSION;
version=str2num([o_ver(1),o_ver(3),o_ver(5)]);
end

if ( ~exist('OCTAVE_VERSION') | version >= 320 )
options=optimset;
options.Display='off';
options.TolFun=1.e-12;
options.MaxFunEvals=10000;
end
for glob_t=tt(2:end)
if ( exist('OCTAVE_VERSION') & version < 320 )
  w = fsolve('beulerfun',glob_y);
else
  w = fsolve(@(w) beulerfun(w),glob_y,options);
end
  u = [u; w.'];
  glob_y = w;
end
t=tt';
clear glob_h glob_t glob_y glob_odefun;
end

function [z]=beulerfun(w)
  global glob_h glob_t glob_y glob_odefun;
  z=w-glob_y-glob_h*glob_odefun(glob_t,w);
end

% Gauss-Newton for nonlinear inverse problem
f=inline('beta1*x/(beta2+x)','beta1','beta2','x');% f(x; beta)
f1=inline('x/(beta2+x)','beta1','beta2','x');% 1st term of Jacobian
f2=inline('-beta1*x/(beta2+x)^2','beta1','beta2','x');% 2nd term of J
% measured data
X = [0.038   0.194  0.425   0.626   1.253   2.500   3.740];
Y = [0.05    0.127  0.094   0.2122  0.2729  0.2665  0.3317];
% initial guess
beta10 = 0.9; beta20 = 0.2;
% setup
m  = length(X);
rk = zeros(m, 1);
Jk = zeros(m, 2);
v  = [beta10, beta20]';
for k=0:5 % iterations
   for i=1:length(X)
     rk(i)    = Y(i) - f (beta10, beta20, X(i));
     Jk(i, 1) =       -f1(beta10, beta20, X(i));
     Jk(i, 2) =       -f2(beta10, beta20, X(i));
   end
  fprintf('%d %0.4g %0.4g %0.4g\n', k, v(1), v(2), norm(rk));
  v = v - (Jk'*Jk) \ (Jk'*rk); % solve G-N equation
  beta10 = v(1);
  beta20 = v(2);
end
% plot the results
h=0.1;
xs = 0; xl = max(X)+0.2;
Xe = xs:h:xl;
Ye = 0*Xe;
for i=1:length(Xe)
   Ye(i) = f(beta10, beta20, Xe(i));
end
plot(X,Y,'ro',Xe,Ye);
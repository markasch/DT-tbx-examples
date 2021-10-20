%--- OMP_binary.m ---%
% Simple OMP example
clear
S = [0,0,1,0,0,1,0,0,1,0]'; % unknown, 3-sparse signal
k = 3; % signal sparsity
% sensing matrix
A = [0  1  1 -1 -1  0 -1  0 -1  0; 
    -1 -1  0  1 -1  0  0 -1  0  1;
     1 -1  1 -1  0 -1  1  1  0  0;
     1  0 -1  0  0  1 -1 -1  1  1;
    -1  0  0  0  1  0  1  0  1 -1;
     0  0 -1 -1 -1  0 -1  1 -1  0];
% measurement
y = A*S; y'
% Step 1:
r = y;
% find argmax (column of X)
for j = 1:10
  vj(j) = r'*A(:,j);
end
maxc1 = find(vj==max(vj))
% residual
gamma = 1;
r = r - gamma*A(:,maxc1);
% Step 2
% find argmax (column of X)
for j = 1:10
  vj(j) = r'*A(:,j);
end
maxc2 = find(vj==max(vj))
r = r - gamma*A(:,maxc2);
% Step 3
% find argmax (column of X)
for j = 1:10
  vj(j) = r'*A(:,j);
end
maxc3 = find(vj==max(vj))
r = r - gamma*A(:,maxc3);
norm(r)
% Recover sparse signal from the argmax columns
Sr = zeros(10,1);
maxc = [maxc1, maxc2, maxc3]
Sr(maxc)=1; Sr'
%
% Better version with a loop
%
r = y;  % initial Residue
O = []; % initialisation of Support Vectors
x = zeros(size(A,2),1); % initialisation of x
for i = 1:k
  c = A'*r;           % Correlation
  [~,ind] = max((c)); % Index with max correlation
  O = [O ind];        % Update support
  Ao = A(:,O');       % Update measurement matrix
  x1 = Ao\y;          % min ||y-Ao*x||
  r = r - gamma*A(:,ind); % ith step residual for 0-1 signal
  %r= y - Ao*x1; % ith step residual for general case
end
x(O')=x1; x'
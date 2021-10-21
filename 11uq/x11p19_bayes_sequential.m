%--- BayesSequential.m ---%
% Sequential Bayesian estimation of 2-parameter system 
% to estimate (x,y) coordinates of a location from noisy measurements
s = [3;5]; % True position, or state of the system.
N = 100;   % Number of observations.
figure(1);clf;
% Create a 100-sample Gaussian noise sequence with a std. dev. of 2.
n = 2*randn(2, 100);
x = zeros(2, 100);
% Add noise to true state to create N observations.
for i = 1:N
  x(:,i) = s + n(:, i);  
end
plot(x(1,:),x(2,:),'k.'); % scatter plot of points	
% Range of possible positions for the system.
Sa = [2:0.05:4];
Sb = [4:0.05:6];
L  = length(Sa);
% Pr is prior distribution
% Po is posterior distribution
Pr = ones(L,L); % Initialize to all ones - uniform prior
Po = ones(L,L);
Pr = Pr/sum(sum(Pr)); % Convert to a pdf by dividing by the sum.
Po = Po/sum(sum(Po)); % Each entry is now 1/L^2.
%Pr=0*Pr;Pr(2,2)=1; % a sharply peaked prior
figure(1); clf; mesh(Po), axis([0 40 0 40 0 0.015]), title('Prior')
[a,b]=find(Po==max(max(Po))); % Find indices where Po has its max.
sest=[Sa(a);Sb(b)];  % The best estimate of the true state to start.
figure(1); clf; figure(2); clf
subplot(211); plot(1,sest(1)); axis([0 N 2 4 ]); hold on;
line([1,N],[s(1),s(1)]); % Draw a line at location of true x-coord.
subplot(212); plot(1,sest(2)); axis([0 N 4 6 ]); hold on;
line([1,N],[s(2),s(2)]); % Draw a line at location of true y-coord.
K=[4,0;0,4]; % covariance matrix for 2-D Gaussian
for n=2:N; % Main loop over observations
  Pr=Po; % Store the posterior as the updated prior.
  m=0*Pr;   
  % Compute likelihood by grid search:
  % - look at each location, assume that the given location is
  %   where object is,
  % - compute likelihood of data x(:,n) assuming 2-d Gaussian noise
  for i = 1:length(Pr)
    for j = 1:length(Pr)
      me = [Sa(i);Sb(j)];
      %Compute likelihood 
      err = x(:,n) - me % observation - model estimation
      m(i,j) = 1/sqrt((2*pi)^2*det(K))*exp(-err'*inv(K)*err/2);           
      m(i,j) = m(i,j) * Pr(i,j); % Bayes: likelihood x prior   
    end;
  end;
  Po = m/sum(sum(m)); % Normalize distribution to make proper pdf.
  figure(1);mesh(Po), axis([0 40 0 40 0 0.015]) % Plot it.
  figure(2);
  [a,b] = find(Po==max(max(Po))); % Max  is most likely location.
  sest  = [Sa(a);Sb(b)];  % Best estimate of the location
  subplot(211);plot(n,sest(1),'k.');axis([0 N 2 4 ])
  subplot(212); plot(n,sest(2),'k.');axis([0 N 4 6 ])
  disp('Key to continue'); pause
end;  
subplot(211); hold off;
subplot(212); hold off;

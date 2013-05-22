% clear, clc,

% set simulation metadata
T       = 1000; % # of time steps
V.dt    = 1/8;  % time step size

% initialize params
P.a     = 1;    % observation scale
P.b     = 0;    % observation bias
tau     = 1.5;    % decay time constant
P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
P.lam   = 0.1;  % firing rate = lam*dt
P.sig   = 0.1;  % standard deviation of observation noise 

% simulate data
N = poissrnd(P.lam*V.dt*ones(T,1)); % simulate spike train
C = filter(1,[1 -P.gam],N);         % calcium concentration
F = P.a*C+P.b + P.sig*randn(T,1);   % observations

% infer spikes
[Nhat Phat] = fast_oopsi(F,V,P);

% plot results
figure(1), clf
tvec=0:V.dt:(T-1)*V.dt;
h(1)=subplot(311); plot(tvec,F); axis('tight'), ylabel('F (au)')
h(2)=subplot(312); plot(tvec,C); axis('tight'), ylabel('C (au)')
h(3)=subplot(313); stem(tvec,N); axis('tight'), ylabel('N (# spikes)')
hold on, plot(tvec,Nhat,'r'), xlabel('time (sec)') % estimated spike train
linkaxes(h,'x'), legend('spike train', 'estimate')
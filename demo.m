clear, clc,

% set simulation metadata
T               = 1500;     % # of time steps
V.dt            = 1/200;    % time step size

% initialize params
P.a     = 200;              % scale
P.b     = 1/P.a;            % bias
tau     = 0.5;              % decay time constant for each cell
P.gam   = 1-V.dt/tau;       % gam
P.lam   = 2.0;              % rate
P.sig   = 50;              % noise

% simulate data
N = poissrnd(P.lam*V.dt*ones(T,1)); % simulate spike train
C = filter(1,[1 -P.gam],N);         % calcium concentration
F = P.a*C+P.b + P.sig*randn(T,1);

% infer spikes
[Nhat Phat] = fast_oopsi(F,V,P);

% plot results
figure(1), clf
h(1)=subplot(311); plot(F); axis('tight')
h(2)=subplot(312); plot(C);
h(3)=subplot(313); stem(N);
hold on, plot(Nhat,'r')
linkaxes(h,'x')
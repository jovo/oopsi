function S = smc_oopsi_forward(F,V,P)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that V.Nspikehist=0 (M is the # of spike history terms).
%
% The backward sampler has standard variance, as approximated typically.
% Each particle has the SAME backwards sampler, initialized at E[h_t]
% this function only does spike history stuff if necessary
%
% the model is F_t = f(C) = alpha C^n/(C^n + k_d) + beta + e_t,
% where e_t ~ N[0, gamma*f(C)+zeta]
%
% Inputs---
% F: Fluorescence
% V: Variables for algorithm to run
% P: initial Parameter estimates
%
% Outputs---
% S: simulation states

%% allocate memory and initialize stuff

fprintf('\nT = %g steps',V.T)
fprintf('\nforward step:        ')
P.kx        = P.k'*V.x;
A.nsmat = ones(V.Nparticles,1) * (0:V.maxspikes);
A.repmat_msp1 = ones(V.maxspikes + 1,1);
IsBernoulli = strcmpi(V.spikegen.name, 'bernoulli');

% extize particle info
S.p     = zeros(V.Nparticles,V.T); %intialize the parameter describing the spiking distribution. in the case of Bernoulli and Poisson, this will also be the rate
S.n     = zeros(V.Nparticles,V.T);                  % extize spike counts
S.next_n = zeros(V.Nparticles,1);
S.C     = P.C_init*ones(V.Nparticles,V.T);          % extize calcium
if V.CaBaselineDrift
    S.Cr     = P.C_init*ones(V.Nparticles,V.T);     % extize resting calcium
    A.epsilon_cr = randn(V.Nparticles,V.T);         % generate noise on Cr
end
S.w_f   = 1/V.Nparticles*ones(V.Nparticles,V.T);    % extize forward weights
S.w_b   = 1/V.Nparticles*ones(V.Nparticles,V.T);    % extize forward weights
S.Neff  = 1/V.Nparticles*ones(1,V.T_o); 			% extize N_{eff}

% preprocess stuff for stratified resampling
ints        = linspace(0,1,V.Nparticles+1);         % generate intervals
diffs       = ints(2)-ints(1);              		% generate interval size
A.U_resamp  = repmat(ints(1:end-1),V.T_o,1)+diffs*rand(V.T_o,V.Nparticles); % resampling matrix

% sample random variables
A.U_sampl   = rand(V.Nparticles,V.T);               % random samples
A.epsilon_c = randn(V.Nparticles,V.T);              % generate noise on C

A.constant_spikedistrib = 0;
if V.Nspikehist>0                                   % if spike histories
    S.h = zeros(V.Nparticles,V.T,V.Nspikehist);     % extize spike history terms
    A.epsilon_h = zeros(V.Nparticles, V.T, V.Nspikehist); % generate noise on h
    for m=1:V.Nspikehist                            % add noise to each h
        A.epsilon_h(:,:,m) = sqrt(P.sig2_h(m))*randn(V.Nparticles,V.T);
    end
else                                                % if not, comput P[n_t] for all t
    S.p = repmat(V.spikegen.nonlinearity(P.kx, P, V)',1,V.Nparticles)';
    if V.freq == 1 && all(S.p(:) == S.p(1))
        A.constant_spikedistrib = 1;
        V.common_ln_n = V.spikegen.logdensity(S.p(:,1)', P, V)'; % compute log(P[n spikes], n = 0,1,2 ... each row is a particle, each column is a spike count
    end
end

% extize stuff needed for conditional sampling
A.n_sampl   = rand(V.Nparticles,V.T);       % generate random number to use for sampling n_t
A.C_sampl   = rand(V.Nparticles,V.T);       % generate random number to use for sampling C_t
A.oney      = ones(V.Nparticles,1);         % vector of ones to call for various purposes to speed things up
A.zeroy     = zeros(V.Nparticles,1);        % vector of zeros

% extize stuff needed for REAL backwards sampling
O.p_o       = zeros(2^(V.freq-1),V.freq);   % extize backwards prob
mu_o      = zeros(2^(V.freq-1),V.freq);   % extize backwards mean
sig2_o    = zeros(1,V.freq);              % extize backwards variance

% extize stuff needed for APPROX backwards sampling
O.p         = zeros(V.freq,V.freq);         % extize backwards prob
O.mu        = zeros(V.freq,V.freq);         % extize backwards mean
O.sig2      = zeros(V.freq,V.freq);         % extize backwards var

% initialize backwards distributions
s           = V.freq;                       % initialize time of next observation
O.p_o(1,s)  = 1;                            % initialize P[F_s | C_s]

[mu_o(1,s) sig2_o(s)]  = init_lik(P,F(s));

O.p(1,s)    = 1;                            % initialize P[F_s | C_s]
O.mu(1,s)   = mu_o(1,s);                  % initialize mean of P[O_s | C_s]
O.sig2(1,s) = sig2_o(s);                  % initialize var of P[O_s | C_s]

if V.freq>1                                 % if intermitent sampling
    for tt=s:-1:2                           % generate spike binary matrix
        A.spikemat(:,tt-1) = repmat([repmat(0,1,2^(s-tt)) repmat(1,1,2^(s-tt))],1,2^(tt-2))';
    end
    nspikes = sum(A.spikemat')';            % count number of spikes at each time step

    for n=0:V.freq-1
        A.ninds{n+1}= find(nspikes==n);     % get index for each number of spikes
        A.lenn(n+1) = length(A.ninds{n+1}); % find how many spikes
    end
end

O = update_moments(V,F,P,S,O,A,s);          % recurse back to get P[O_s | C_s] before the first observation

%% do the particle filter

for t=V.freq+1:V.T-V.freq
    if IsBernoulli && V.freq > 1 && ~V.CaBaselineDrift && V.condsamp
        S = cond_sampler_intermittent_bernoulli_nodrift(V,F,P,S,O,A,t,s);
    elseif V.condsamp==0 || (F(s+V.freq)-P.beta)/P.alpha>0.98 || (F(s+V.freq)-P.beta)/P.alpha<0 % prior sampler with drift when gaussian approximation is not good
        S = prior_sampler(V,F,P,S,A,t);     % use prior sampler when gaussian approximation is not good
    else % otherwise use conditional sampler
        S = cond_sampler_nonintermittent(V,F,P,S,O,A,t,s);
    end

    S.C(:,t)=S.next_C;                      % done to speed things up for older
    if V.CaBaselineDrift
        S.Cr(:,t) = S.next_Cr;
    end
    S.n(:,t)=S.next_n;                      % matlab, having issues with from-function
    if(isfield(S,'next_w_f'))               % update calls to large structures
        S.w_f(:,t)=S.next_w_f;
    else
        S.w_f(:,t)=S.w_f(:,t-1);
    end

    if V.Nspikehist>0                                % update S.h & S.p
        for m=1:V.Nspikehist, S.h(:,t,m)=S.h_new(:,1,m); end
        S.p(:,t)=S.p_new;
    end

    % at observations
    if mod(t,V.freq)==0
        % STRAT_RESAMPLE -- THIS WAS CAUSING TROUBLE IN OLDER MATLAB WITH
        % MANAGING LARGE ARRAYS IN S-STRUCTURE WITHIN THE FUNCTION CALL
        % S = strat_resample(V,S,t,A.U_resamp); % stratified resample

        Nresamp=t/V.freq;                                       % increase sample counter
        S.Neff(Nresamp)  = 1/sum(S.w_f(:,t).^2);                % store N_{eff}

        % if weights are degenerate or we are doing prior sampling then resample
        if S.Neff(Nresamp) < V.Nparticles/2 || V.condsamp==0
            [foo,ind]   = histc(A.U_resamp(Nresamp,:),[0  cumsum(S.w_f(:,t))']);
            [ri,ri]     = sort(rand(V.Nparticles,1));                    % these 3 lines stratified resample
            ind         = ind(ri);

            S.p(:,t-V.freq+1:t)   = S.p(ind,t-V.freq+1:t);      % resample probabilities (necessary?) (yes, i think we need these for the backward step -- dg)
            S.n(:,t-V.freq+1:t)   = S.n(ind,t-V.freq+1:t);      % resample spikes
            S.C(:,t-V.freq+1:t)   = S.C(ind,t-V.freq+1:t);      % resample calcium
            if V.CaBaselineDrift
                S.Cr(:,t-V.freq+1:t)   = S.Cr(ind,t-V.freq+1:t);      % resample resting calcium
            end
            S.w_f(:,t-V.freq+1:t) = 1/V.Nparticles*ones(V.Nparticles,V.freq); % reset weights
            if V.Nspikehist>0                                   % if spike history terms
                S.h(:,t-V.freq+1:t,:) = S.h(ind,t-V.freq+1:t,:);% resample all h's
            end
        end %function

        O = update_moments(V,F,P,S,O,A,t);                      % estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
        s = t;                                                  % store time of last observation

    end

    if mod(t,100)==0
        if t<1000                                               % print # of observations
            fprintf('\b\b\b%d',t)
        elseif t<10000
            fprintf('\b\b\b\b%d',t)
        elseif t<100000
            fprintf('\b\b\b\b\b%d',t)
        end
    end

end %for time loop

end %conditional sampler

%% initialize likelihood
function [mu1 sig1] = init_lik(P,F)
% get the mean (mu1) and variance (sig1) for P[C_t | F_t]
% compute mean
num = P.k_d * (P.beta - F);
denom = F - P.beta - P.alpha;
ginv = (num ./ denom) .^ (1 / P.n);
mu1 = ginv;

% compute var
if mu1>0 && imag(mu1)==0
    mu_nm1power = mu1 .^ (P.n - 1);
    gprime = P.alpha * P.n * mu_nm1power * P.k_d ./ ((mu1 .* mu_nm1power + P.k_d).^2);

    SCa_at_mu1 = (F - P.beta) / P.alpha;
    sig1 = (P.gamma * SCa_at_mu1 + P.zeta) ./ (gprime^2);
else
    mu1=0;
    sig1=0;
end

end %init_lik

%% update moments

function O = update_moments(V,F,P,S,O,A,t)
%%%% maybe make a better proposal for epi

s           = V.freq;                   % find next observation time

[mu1 sig1]  = init_lik(P,F(t+s));

O.p(1,s)    = 1;                        % initialize P[F_s | C_s]
O.mu(1,s)   = mu1;                      % initialize mean of P[O_s | C_s]
O.sig2(1,s) = sig1;                     % initialize var of P[O_s | C_s]

if s > 1
    if V.CaBaselineDrift
        error('intermittent sampling + baseline drift is not yet implemented');
    end
    if ~strcmpi(V.spikegen.name, 'bernoulli')
        error('intermittent sampling is only implemented for bernoulli spiking');
    end

    mu_o(1,s) = mu1;                      % initialize mean of P[O_s | C_s]
    sig2_o(s) = sig1;                     % initialize var of P[O_s | C_s]


    if V.Nspikehist>0
        hhat        = zeros(V.freq,V.Nspikehist);                   % extize hhat
        phat        = zeros(1,V.freq+1);                            % extize phat

        hs          = S.h(:,t,:);                                   % this is required for matlab to handle a m-by-n-by-p matrix
        h(:,1:V.Nspikehist)= hs(:,1,1:V.Nspikehist);                % this too
        hhat(1,:)   = sum(repmat(S.w_f(:,t),1,V.Nspikehist).*h,1);  % initialize hhat
        phat(1)     = sum(S.w_f(:,t).*S.p(:,t),1);                  % initialize phat
    end

    if V.Nspikehist>0
        for tt=1:s
            % update hhat
            for m=1:V.Nspikehist                                    % for each spike history term
                hhat(tt+1,m)=(1-P.g(m))*hhat(tt,m)+phat(tt);
            end
            y_t         = P.kx(tt+t)+P.omega'*hhat(tt+1,:)';        % input to neuron
            phat(tt+1)  = 1-exp(-exp(y_t)*V.dt);                    % update phat
        end
    else
        phat  = 1-exp(-exp(P.kx(t+1:t+s)')*V.dt);                   % update phat
    end

    for tt=s:-1:2

        O.p_o(1:2^(s-tt+1),tt-1)    = repmat(O.p_o(1:2^(s-tt),tt),2,1).*[(1-phat(tt))*ones(1,2^(s-tt)) phat(tt)*ones(1,2^(s-tt))]';
        mu_o(1:2^(s-tt+1),tt-1)   = (1-P.a)^(-1)*(repmat(mu_o(1:2^(s-tt),tt),2,1)-P.A*A.spikemat(1:2^(s-tt+1),tt-1)-P.a*P.C_0);     %mean of P[O_s | C_k]
        sig2_o(tt-1)              = (1-P.a)^(-2)*(P.sig2_c+sig2_o(tt)); % var of P[O_s | C_k]

        for n=0:s-tt+1
            nind=A.ninds{n+1};
            O.p(n+1,tt-1)   = sum(O.p_o(nind,tt-1));
            ps          = (O.p_o(nind,tt-1)/O.p(n+1,tt-1))';
            O.mu(n+1,tt-1)  = ps*mu_o(nind,tt-1);
            O.sig2(n+1,tt-1)= sig2_o(tt-1) + ps*(mu_o(nind,tt-1)-repmat(O.mu(n+1,tt-1)',A.lenn(n+1),1)).^2;
        end
    end

    if s==2
        O.p     = O.p_o;
        O.mu    = mu_o;
        O.sig2  = repmat(sig2_o,2,1);
        O.sig2(2,2) = 0;
    end

    while any(isnan(O.mu(:)))              % in case ps=0/0, which yields a NaN, approximate mu and sig
        O.mu(1,:)   = mu_o(1,:);
        O.sig2(1,:) = sig2_o(1,:);
        ind         = find(isnan(O.mu));
        O.mu(ind)   = O.mu(ind-1)-P.A;
        O.sig2(ind) = O.sig2(ind-1);
    end
end
O.p=O.p+eps; % such that there are no actual zeros
end %function UpdateMoments

%% particle filtering using the prior sampler
function S = prior_sampler(V,F,P,S,A,t)
S = calc_next_spikedist_params(V,P,S,A,t);
if V.use_true_n
    S.next_n    = V.true_n(:,t);
else
    den = V.spikegen.density(S.p_new', P, V);
    for pind = 1:V.Nparticles
        cden = cumsum([0 den(:,pind)']);
        cden(end) = 1 + eps;
        S.next_n(pind,1) = find(A.U_sampl(pind,t) > cden,1,'last') - 1;
    end
end
if V.CaBaselineDrift
    S.next_Cr = S.Cr(:,t-1) + sqrt(P.sig2_cr) * A.epsilon_cr(:,t);
    S.next_C        = (1-P.a) * S.C(:,t-1) + P.A * S.next_n + P.a * S.Cr(:,t) + sqrt(P.sig2_c) * A.epsilon_c(:,t);% sample C
else
    S.next_C        = (1-P.a) * S.C(:,t-1) + P.A *S.next_n + P.a*P.C_0 + sqrt(P.sig2_c) * A.epsilon_c(:,t);% sample C
end
% get weights at every observation          %THIS NEEDS FIX FOR EPI DATA
if mod(t,V.freq)==0
    S_mu        = Hill_v1(P,S.next_C);
    F_mu        = P.alpha*S_mu+P.beta;      % compute E[F_t]
    F_var       = P.gamma*S_mu+P.zeta;      % compute V[F_t]
    %%%% this must also change for epi
    ln_w        = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2;% compute log of weights
    ln_w        = ln_w-max(ln_w);           % subtract the max to avoid rounding errors
    w           = exp(ln_w);                % exponentiate to get actual weights
    %     error('forgot to include the previous weight in this code!!!!')
    %     break
    S.next_w_f  = w/sum(w);                 % normalize to define a legitimate distribution
end

end

function S = calc_next_spikedist_params(V,P,S,A,t) % if spike histories, sample h and update p
if V.Nspikehist>0
    S.h_new=zeros(size(S.n,1),1,V.Nspikehist);
    for m=1:V.Nspikehist
        S.h_new(:,1,m)=(1-P.g(m))*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end
    % update rate and sample spikes
    hs          = S.h_new;              % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:V.Nspikehist)  = hs(:,1,1:V.Nspikehist);        % this too
    y_t         = P.kx(t)+P.omega'*h';  % input to neuron
    S.p_new     = V.spikegen.nonlinearity(y_t, P, V);
    %S.p_new     = 1-exp(-exp(y_t)*V.dt);% update rate for those particles with y_t<0
    S.p_new     = S.p_new(:);
else
    S.p_new         = S.p(:,t);
end
end

function S = cond_sampler_intermittent_bernoulli_nodrift(V,F,P,S,O,A,t,s)
S = calc_next_spikedist_params(V,P,S,A,t);
k   = V.freq-(t-s)+1;
A.repmat_ktimes_mat = ones(k,1);
m2  = O.mu(1:k,t-s);                        % mean of P[O_s | C_k] for n_k=1 and n_k=0
v2  = O.sig2(1:k,t-s);                      % var of P[O_s | C_k] for n_k=1 and n_k=0
prev_C = S.C(:,t-1);

% compute P[n_k | h_k]
ln_n    = [log(S.p_new) log(1-S.p_new)];    % compute [log(spike) log(no spike)]

m0  = (1-P.a) * prev_C + P.a *  P.C_0;         % mean of P[C_k | C_{t-1}, n_k=0]
m1  = (1-P.a) * prev_C + P.a *  P.C_0 + P.A ;     % mean of P[C_k | C_{t-1}, n_k=1]

v = P.sig2_c + v2;
v = v(:,A.oney)';                            % var of G_n(n_k | O_s) for n_k=1 and n_k=0

%is the -0.5*log(2*pi.*v) term necessary ???
ln_G0 = -0.5*log(2*pi * v) - 0.5 * (m0(:,A.repmat_ktimes_mat) - m2(:,A.oney)').^2 ./ v;
ln_G1 = -0.5*log(2*pi * v) - 0.5 * (m1(:,A.repmat_ktimes_mat) - m2(:,A.oney)').^2 ./ v;

mx = max(max(ln_G0,[],2),max(ln_G1,[],2));% get max of these
mx = mx(:,A.repmat_ktimes_mat);

%why calculate exp when we're soon going to multiply by something and then take logs again ???
G1  = exp(ln_G1-mx);                        % norm dist'n for n=1;
M1  = G1*O.p(1:k,t-s);                      % times prob of n=1

G0  = exp(ln_G0-mx);                        % norm dist'n for n=0;
M0  = G0*O.p(1:k,t-s);                      % times prob n=0

ln_G    = [log(M1) log(M0)];                % ok, now we actually have the gaussians

% compute q(n_k | h_k, O_s)
ln_q_n  = ln_n + ln_G;                      % log of sampling dist
mx      = max(ln_q_n,[],2);                 % find max of each row
mx2     = [mx mx];                          % matricize
q_n     = exp(ln_q_n-mx2);                  % subtract max to ensure that for each row, there is at least one positive probability, and exponentiate
sq_n = sum(q_n,2);
q_n     = q_n ./ [sq_n sq_n];               % normalize to make this a true sampling distribution (ie, sums to 1)

% sample n
S.next_n= A.n_sampl(:,t)<q_n(:,1);          % sample n
sp      = S.next_n==1;                      % store index of which samples spiked
nosp    = S.next_n==0;                      % and which did not

% sample C; intermittent, sample from mixture
% first sample component
if(isempty(find(sp,1))), sp_i=[];       % handler for empty spike trains
else [fo,sp_i]   = histc(A.C_sampl(sp,t),[0  cumsum(O.p(1:k-1,t-s))'/sum(O.p(1:k-1,t-s))]); end
if(isempty(find(nosp,1))), nosp_i=[];   % handle for saturated spike trains
else [fo,nosp_i] = histc(A.C_sampl(nosp,t),[0  cumsum(O.p(1:k,t-s))'/sum(O.p(1:k,t-s))]); end

v       = O.sig2(1:k,t-s);              % get var of each component
v(sp)   = v(sp_i);                      % if particle spiked, then use variance of spiking
v(nosp) = v(nosp_i);                    % o.w., use non-spiking variance

m       = O.mu(1:k,t-s);                % get mean of each component
m(sp)   = m(sp_i);                      % if particle spiked, use spiking mean
m(nosp) = m(nosp_i);                    % o.w., use non-spiking mean

v_c         = (1./v + 1/P.sig2_c).^(-1);      %variance of dist'n for sampling C
m_givenspiking = (1-P.a) * S.C(:,t-1) + P.A * S.next_n + P.a * C_baseline; %expected calcium given spiking and previous calcium value
m_c         = v_c .* (m ./ v + m_givenspiking / P.sig2_c);%mean of dist'n for sampling C
S.next_C    = normrnd(m_c, sqrt(v_c));       % sample C

% update weights
if mod(t,V.freq)==0                         % at observations compute P(O|H)
    S_mu        = Hill_v1(P,S.next_C);
    if V.scan==0                            % when doing epi, also add previous time steps
        for tt=s+1:t-1, S_mu=S_mu+Hill_v1(P,S.C(s+tt)); end
    end
    F_mu        = P.alpha*S_mu+P.beta;      % compute E[F_t]
    F_var       = P.gamma*S_mu+P.zeta;      % compute V[F_t]
    %%%% log_PO_H must change for epi
    log_PO_H    = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2; % compute log of P[F | H]
else                                        % when no observations are present, P[F |H] is equal for all particles
    log_PO_H    = (1/V.Nparticles)*A.oney;
end

log_n           = A.oney;                   % extize log sampling spikes
log_n(sp)       = log(S.p_new(sp));         % compute log P(spike)
log_n(nosp)     = log(1-S.p_new(nosp));     % compute log P(no spike)
log_C_Cn        = -0.5*(S.next_C-((1-P.a)*S.C(:,t-1)+P.A*S.next_n+P.a* C_baseline)).^2/P.sig2_c;%log P[C_k | C_{t-1}, n_k]

log_q_n         = A.oney;                   % initialize log q_n
log_q_n(sp)     = log(q_n(sp,1));           % compute what was the log prob of sampling a spike
log_q_n(nosp)   = log(1-q_n(nosp,1));       % or sampling no spike
log_q_C         = -0.5*(S.next_C - m_c).^2./v_c;% log prob of sampling the C_k that was sampled

log_quotient    = log_PO_H + log_n + log_C_Cn + log_q_n - log_q_C; %note that terms from baseline drift cancel
sum_logs        = log_quotient+log(S.w_f(:,t-1));   % update log(weights)
w               = exp(sum_logs-max(sum_logs));      % exponentiate log(weights)
S.next_w_f      = w./sum(w);                        % normalize such that they sum to unity
if any(isnan(w)) %, Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
    warning('smc:weights','some weights are NaN')
    keyboard,
end

end

function S = cond_sampler_nonintermittent(V,F,P,S,O,A,t,s)
%4 steps: 1)sample n, 2)sample Cr (if drift allowed), 3)sample C, 4)update weights
mCF  = O.mu;                        % mean of P[F | C_k]
vCF  = O.sig2;                      % var of P[F | C_k]
%       ----        STEP 1: Sample n       ----
if A.constant_spikedistrib %the parameters of the spike count distributions are constant for all particles and time steps
    ln_n = V.common_ln_n;  % log(P[n spikes], n = 0,1,2 ... each row is a particle, each column is a spike count
else %the parameters of the spike count distributions can vary over time and particles, due to a stimulus or spike history terms
    S = calc_next_spikedist_params(V,P,S,A,t);
    ln_n = V.spikegen.logdensity(S.p_new', P, V)'; % compute log(P[n spikes], n = 0,1,2 ... each row is a particle, each column is a spike count
end
prev_C = S.C(:,t-1);
if V.CaBaselineDrift
    prev_Cr = S.Cr(:,t-1);
    mFCrn = (mCF - prev_C(:,A.repmat_msp1) * (1 - P.a) - P.A * A.nsmat) / P.a;  %mFCrn(i,j) is mean P[F | Cr_t, n_t = j-1] for particle i. derived by integrating out C using the rule for a gaussian product integrand
    vFCr = (vCF + P.sig2_c) / (P.a^2);                                                      %scalar. var of P[F | Cr_t]
    v = vFCr + P.sig2_cr;                                                      %variance of the gaussian in the formula describing P[F | n] i.e. G_n(n_k | F)
    ln_G = - 0.5 * (mFCrn - prev_Cr(:,A.repmat_msp1)).^2 / v;    %common terms ignored: - 0.5*log(2*pi*v) - log(P.a)
else
    mCC = (1-P.a) * prev_C(:,A.repmat_msp1) + P.a * P_C0 + P.A * A.nsmat;       %mCC(i,j) mean of P[C_k | C_{t-1}, n_k = j-1] for particle i
    v = P.sig2_c + vCF;
    ln_G = -0.5*log(2*pi * v) - 0.5 * (mCC - mCF).^2 ./ v; %log G_n(n_k | O_s) for each possible value of n_k:
end
ln_q_n  = ln_n + ln_G;                      % log of sampling dist q(n_k | h_k, O_s)
mx      = max(ln_q_n,[],2);                 % find max of each row
q_n     = exp(ln_q_n - mx(:,A.repmat_msp1));                  % subtract max to ensure that for each row, there is at least one positive probability, and exponentiate
sq_n = sum(q_n,2);
q_n     = q_n ./ sq_n(:,A.repmat_msp1);               % normalize to make this a true sampling distribution (ie, sums to 1)
cden_mat = cumsum([A.zeroy q_n(:,1:end-1)], 2);
S.next_n = sum(A.U_sampl(:, t * A.repmat_msp1) > cden_mat,2) - 1;
%       ----        STEP 2: Sample Cr if baseline is allowed to drift       ----
if V.CaBaselineDrift %sample Cr
    mFCrn = (mCF - prev_C * (1 - P.a) - P.A * S.next_n) / P.a; %mean of P[F | Cr_t, n_t = S.next_n]
    v_cr = (1 ./ vFCr + 1 / P.sig2_cr).^(-1);   % variance of distribution for sampling Cr
    m_cr = v_cr .* (prev_Cr / P.sig2_cr + mFCrn ./ vFCr); %mean of distribution for sampling Cr
    S.next_Cr = m_cr + sqrt(v_cr) * A.epsilon_cr(:,t);
    C_baseline = S.next_Cr;
else
    C_baseline = P.C_0;
end
%       ----        STEP 3: Sample C       ----
v_c         = 1/(1/vCF + 1/P.sig2_c);      %variance of dist'n for sampling C
m_givenspiking = (1-P.a) * S.C(:,t-1) + P.A * S.next_n + P.a * C_baseline; %expected calcium given spiking and previous calcium value
m_c         = v_c * (mCF / vCF + m_givenspiking / P.sig2_c);%mean of dist'n for sampling C
S.next_C = m_c + sqrt(v_c) * A.epsilon_c(:,t);
%       ----        STEP 4: update weights      ----
S_mu        = Hill_v1(P,S.next_C);
F_mu        = P.alpha*S_mu+P.beta;      % compute E[F_t | C_t]
F_var       = P.gamma*S_mu+P.zeta;      % compute V[F_t | C_t]
log_PF_H    = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2; % compute log of P[F | H] !!!MUST CHANGE FOR EPI!!!

direct_ind = (1:V.Nparticles)' + S.next_n * V.Nparticles; %index to access ln_n(pind, S.next_n(pind + 1) for pind = 1:V.nparticles
log_n = ln_n(direct_ind);
log_C_Cn        = -0.5*(S.next_C-((1-P.a)*S.C(:,t-1)+P.A*S.next_n+P.a* C_baseline)).^2/P.sig2_c;%log P[C_k | C_{t-1}, n_k]

log_q_n = log(q_n(direct_ind));             % compute what was the log prob of sampling the number of spikes we sampled, for each particle.
log_q_C         = -0.5*(S.next_C - m_c).^2./v_c;% log prob of sampling the C_k that was sampled
if V.CaBaselineDrift
    log_Cr_Cr = -0.5*(S.next_Cr - S.Cr(:,t-1)).^2 / P.sig2_cr;
    log_q_Cr =  -0.5*(S.next_Cr - m_cr).^2./v_cr ;% log prob of sampling the Cr_k that was sampled
else
    log_Cr_Cr = 0;
    log_q_Cr = 0;
end
log_quotient    = log_PF_H + log_n + log_C_Cn + log_Cr_Cr - log_q_n - log_q_C - log_q_Cr; %note that terms from baseline drift cancel
sum_logs        = log_quotient+log(S.w_f(:,t-1));   % update log(weights)
w               = exp(sum_logs-max(sum_logs));      % exponentiate log(weights)
S.next_w_f      = w./sum(w);                        % normalize such that they sum to unity
if any(isnan(w)) %, Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
    warning('smc:weights','some weights are NaN')
    keyboard;
end
end
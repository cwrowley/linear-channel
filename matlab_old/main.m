% main - driver for linear channel code
%
% compute balanced truncation, BPOD, DMD, and ERA
% evaluate error norms
%
% Modified so that state is always in (velocity, vorticity) form
% and inner product does not have the Laplacian in it
% Many of the dependent files will be changed in this directory
clear all
close all

Nz = 16;
Ny = 32;
Re = 1000;

nmodes = 30;    % number of modes to keep
ntrunc = 10;    % number of modes to keep in reduced-order model
projrank = 10;  % rank of output projection in approximate BT
dt = 0.05;
dt2 = 0.05;     % smaller timestep, for computing dmd to avoid d2c errors % keep them the same for now
Tf = 1000;%1000; %1000
tp = 0 : dt : Tf;
tp2 = 0:dt2:Tf;

bal_trun = 1;   % 1- compute exact balanced truncation, 0- don't
inf_norm_flag = 1;  % 1- compute infinity norms; 0: don't
icFlag = 1; % do optimal growth initial condition. Otherwise, as in original code

if ~exist('chan','var')
    %Aadj is now defined in consistant coords
    [Amat, Aadj, Bmat, Mmat, Npts, y, z] = define_eqns(Ny, Nz, Re);
    Cmat = eye(2*Npts);
    %%
    if icFlag
        Tgrowths =32.9;% optimal %[32.9:0.01:33.1];%[1 30 50 200];%[0.1:0.1:100]; %optimal time
        growthEnergies = zeros(2*Npts,length(Tgrowths));
        growthVecs = zeros(2*Npts,length(Tgrowths));
        
        
        for kk = 1:length(Tgrowths)
            %k = 200;
            expAmat = expm(Amat*Tgrowths(kk));
            [U,S,V] = svd(expAmat);
            % svd using randomized projections if state is too large
            %[U,S,V] =rsvdloc(expAmat,k);
            growthEnergies(:,kk) = diag(S); %note that this needs to be squared to be have energy units
            growthVecs(:,kk) = V(:,1);
        end
        
        [~,growthInd] = max(growthEnergies(1,:));
        
        Bmat = growthVecs(:,growthInd);
        % Bmat = sum(growthVecs,2);
    end
    
    %%
    chan = ss(Amat,Bmat,Cmat,0);
    
    %    if (bal_trun ==1)
    disp('Computing minimal realization')
    tic
    chanmin = minreal(chan); % use for computing norms
    t = toc; disp(sprintf('    time = %.2f',t));
    %    end
end

%%
%% Compute Gramians
%%
if ~exist('Wc','var')
    disp('Computing controllability Gramian...')
    tic
    Wc = gram(chan, 'c');
    t = toc; disp(sprintf('    time = %.2f',t));
    
    disp('Computing observability Gramian...')
    tic
    Wo = gram(chan, 'o');
    t = toc; disp(sprintf('    time = %.2f',t));
    
    disp('Computing SVD of Gramians...')
    tic
    [Uc, sigc] = svd(Wc); sigc = diag(sigc);
    [Uo, sigo] = svd(Wo); sigo = diag(sigo);
    t = toc; disp(sprintf('    time = %.2f',t));
    
    % single-output model: most controllable mode (first POD mode)
    Cmat1 = Uc(:,1)';
    chan1 = ss(Amat, Bmat, Cmat1, 0);
end

% Gramians for single-output system
if ~exist('Wo1','var')
    % controllability Gramian is the same
    Wc1 = Wc;
    
    disp('Computing observability Gramian for SISO system...')
    tic
    Wo1 = gram(chan1, 'o');
    t = toc; disp(sprintf('    time = %.2f',t));
    
    disp('Computing SVD of Gramians...')
    tic
    [Uo1, sigo1] = svd(Wo1); sigo1 = diag(sigo1);
    t = toc; disp(sprintf('    time = %.2f\n',t));
end

%%
%% Exact balanced truncation
%%
if ~exist('T','var')
    
    if (bal_trun == 1)
        disp('Computing exact balanced truncation, single output')
        tic
        [T1, sig1, Tinv1] = sysbalcwr(Amat, Bmat, Cmat1);
        t = toc; disp(sprintf('    time = %.2f\n',t));
        
        disp('Computing exact balanced truncation, full state output')
        tic
        [T, sig, Tinv] = sysbalcwr(Amat, Bmat, Cmat);
        t = toc;
        disp(sprintf('   time = %.2f\n',t));
    end
end

%%
%% POD modes
%%
if ~exist('T_pod','var')
    disp('Computing primal snapshots and output projection')
    tic
    % NOTE: following routine modified NOT to use Mmat as inner prod
    % SD: Now removed Mmat entirely.
    [T_pod, sig_pod, primal] = compute_pod(chan,tp);
    figure
    semilogy(sig_pod,'rx-')
    % output projection matrix
    Phi = T_pod(:,1:projrank);
    P = Phi*Phi';
    t = toc; disp(sprintf('    time = %.2f\n',t));
end

%%
%% Adjoint data
%%
if ~exist('dualdat','var')
    disp('Running adjoint simulations')
    tic
    dualdat = compute_adjoint_data(Aadj,Phi,tp);
    t = toc; disp(sprintf('    time = %.2f\n',t));
end
if ~exist('T_snap10','var')
    nt = length(tp);
    disp('Computing approximate balancing transformation (SISO)')
    tic
    [T1_snap, sig1_snap, Tinv1_snap] = snapshotbal4(primal, dualdat(:,1:nt), tp, nmodes);
    t = toc; disp(sprintf('    time = %.2f\n',t));
    
    disp('Computing approximate balancing transformation (5-mode projection)')
    tic
    [T_snap5, sig_snap5, Tinv_snap5] = snapshotbal4(primal, dualdat(:,1:5*nt), tp, nmodes);
    t = toc; disp(sprintf('    time = %.2f\n',t));
    
    disp('Computing approximate balancing transformation (10-mode projection)')
    tic
    [T_snap10, sig_snap10, Tinv_snap10] = snapshotbal4(primal, dualdat, tp, nmodes);
    t = toc; disp(sprintf('    time = %.2f\n',t));
end


%%
%% Compute error norms
%%
%if ~exist('podnorm_2','var')
if (bal_trun == 1)
    disp('Computing 2-norm of full model')
    tic, channorm_2 = norm(chanmin, 2); t = toc;
    disp(sprintf('     norm: %.3f  (%3.2f sec)', channorm_2, t))
    
    disp('Computing inf-norm of full model')
    tic, channorm_inf = norm(chanmin, inf); t = toc;
    disp(sprintf('     norm: %.3f  (%3.2f sec)', channorm_inf, t))
end
%% study growth of rank
%{
for ii = 1:length(tp)-1;
    Rank(ii) = rank(primal(:,2:(ii+1)));
end
figure
plot(Rank)
%}


%% Inital computation of dmd model with a large initial truncation
rDMDinit = 15; %for much above 30, system is slightly unstable
%try running at smaller timestep to mitigate time-discretization errors
[T_pod2, sig_pod2, primal2] = compute_pod(chan,tp2);
PODcoeffs = T_pod2'*primal2;
[Phi, Lambda, U, S, V, Atilde] =  DMDext(primal2(:,1:end),rDMDinit);
chan_dmd = ss(Atilde,U'*Bmat*dt2,U,0,dt2);
chan_dmdcts = d2c(chan_dmd);
Amatdmd = chan_dmdcts.a;
Bmatdmd = chan_dmdcts.b;
Cmatdmd = chan_dmdcts.c;

[Tdmd, sigdmd, Tinvdmd] = sysbalcwr(Amatdmd, Bmatdmd, Cmatdmd);

%% do era with same parameters as dmd, (this gave identical results, even for high values of m).

YY = permute( T_pod(:,1:rDMDinit)'*primal(:,1:end),[1,3,2]); % output projection
%(mr,n) are size of Hankel Matrix for ERA, choose to be such that
m = 5;
n = length(tp)-m;
[Aera,Bera,Cera,Dera,HSVs] = era(YY,m,n,1,rDMDinit,rDMDinit); %nout = r
chan_era = ss(Aera,Bera*dt,T_pod(:,1:rDMDinit)*Cera,0,dt); % *dt due to convention of dss system
chan_eracts = d2c(chan_era);
Amatera = chan_eracts.a;
Bmatera = chan_eracts.b;
Cmatera = chan_eracts.c;

[Tera, sigera, Tinvera] = sysbalcwr(Amatera, Bmatera, Cmatera);


%%
rMax =10;
for r=1:rMax
    %%
    %% reduced order models
    %%
    disp(sprintf('Computing reduced order models, dimension %d',r))
    tic
    chan_bt     = reduced_model(Amat, Bmat, Cmat, T,        Tinv, r);
    chan_btdmd  = reduced_model(Amatdmd, Bmatdmd, Cmatdmd, Tdmd,        Tinvdmd, r);
    chan_btera  = reduced_model(Amatera, Bmatera, Cmatera, Tera,        Tinvera, r);
    chan_pod    = reduced_model(Amat, Bmat, Cmat, T_pod,    T_pod', r);
    %use adjoint modes to project pod modes. Doesn't seem to work
    chan_podmod = reduced_modelmod(Amat, Bmat, Cmat, T_pod,        Tinv_snap10, r);
    
    models(r).chan_bt =chan_bt;
    models(r).chan_btdmd =chan_btdmd;
    models(r).chan_btera =chan_btera;
    models(r).chan_pod =chan_pod;
    models(r).chan_podmod =chan_podmod;
    
    chan_snap5  = reduced_model(Amat, Bmat, Cmat, T_snap5,  Tinv_snap5, r);
    chan_snap10 = reduced_model(Amat, Bmat, Cmat, T_snap10, Tinv_snap10, r);
    models(r).chan_snap5 =chan_snap5;
    models(r).chan_snap10 =chan_snap10;
    t = toc; disp(sprintf('    time = %.2f\n',t));
    
    % do truncated DMD
    [Phi, Lambda, U, S, V, Atilde] =  DMDext(primal(:,1:end),r);
    chan_dmd = ss(Atilde,U'*Bmat*dt,U,0,dt);
    chan_dmdcts = d2c(chan_dmd);
    models(r).chan_dmd = chan_dmdcts;
    
    YY = permute( T_pod(:,1:r)'*primal(:,1:end),[1,3,2]); % output projection
    %(mr,n) are size of Hankel Matrix for ERA, choose to be such that
    m = 50;
    n = length(tp)-m;
    [Aera,Bera,Cera,Dera,HSVs] = era(YY,m,n,1,r,r); %nout = r
    chan_era = ss(Aera,Bera*dt,T_pod(:,1:r)*Cera,0,dt); % *dt due to convention of dss system
    chan_eracts = d2c(chan_era);
    models(r).chan_era = chan_eracts;
    
    % do ERA with a fixed 10 modes
    
    YY2 = permute( T_pod(:,1:10)'*primal(:,1:end),[1,3,2]); % output projection
    %(mr,n) are size of Hankel Matrix for ERA, choose to be such that
    m = 50;
    n = length(tp)-m;
    [Aera,Bera,Cera,Dera,HSVs] = era(YY2,m,n,1,10,r); %nout = r
    chan_era10 = ss(Aera,Bera*dt,T_pod(:,1:10)*Cera,0,dt); % *dt due to convention of dss system
    chan_eracts10 = d2c(chan_era10);
    models(r).chan_era10 = chan_eracts10;
    
    
    disp('Computing 2-norms')
    tic, btnorm_2(r) = norm(chanmin - chan_bt,2); t=toc;
    disp(sprintf('     BT: %.3f  (%3.2f sec)', btnorm_2(r), t))
    tic, podnorm_2(r) = norm(chanmin - chan_pod,2); t=toc;
    disp(sprintf('    pod: %.3f  (%3.2f sec)', podnorm_2(r), t))
    tic, snap5norm_2(r) = norm(chanmin - chan_snap5,2); t=toc;
    disp(sprintf('  snap5: %.3f  (%3.2f sec)', snap5norm_2(r), t))
    tic, snap10norm_2(r) = norm(chanmin - chan_snap10,2); t=toc;
    disp(sprintf(' snap10: %.3f  (%3.2f sec)\n', snap10norm_2(r), t))
    tic, dmdnorm_2(r) = norm(chanmin - chan_dmdcts,2); t=toc;
    disp(sprintf(' dmd: %.3f  (%3.2f sec)\n', dmdnorm_2(r), t))
    tic, dmdbtnorm_2(r) = norm(chanmin - chan_btdmd,2); t=toc;
    disp(sprintf(' dmdbt: %.3f  (%3.2f sec)\n', dmdbtnorm_2(r), t))
    tic, erabtnorm_2(r) = norm(chanmin - chan_btera,2); t=toc;
    disp(sprintf(' erabt: %.3f  (%3.2f sec)\n', erabtnorm_2(r), t))
    tic, eranorm_2(r) = norm(chanmin - chan_eracts,2); t=toc;
    disp(sprintf(' era: %.3f  (%3.2f sec)\n', eranorm_2(r), t))
    tic, era10norm_2(r) = norm(chanmin - chan_eracts10,2); t=toc;
    disp(sprintf(' era10mode: %.3f  (%3.2f sec)\n', era10norm_2(r), t))
    tic, podmodnorm_2(r) = norm(chanmin - chan_podmod,2); t=toc;
    disp(sprintf(' podmod: %.3f  (%3.2f sec)\n', podmodnorm_2(r), t))
    if inf_norm_flag
        disp('Computing infinity norms')
        tic, btnorm_inf(r) = norm(chanmin - chan_bt,inf); t=toc;
        disp(sprintf('     BT: %.3f  (%3.2f sec)', btnorm_inf(r), t))
        tic, podnorm_inf(r) = norm(chanmin - chan_pod,inf); t=toc;
        disp(sprintf('    pod: %.3f  (%3.2f sec)', podnorm_inf(r), t))
        tic, snap5norm_inf(r) = norm(chanmin - chan_snap5,inf); t=toc;
        disp(sprintf('  snap5: %.3f  (%3.2f sec)', snap5norm_inf(r), t))
        tic, snap10norm_inf(r) = norm(chanmin - chan_snap10,inf); t=toc;
        disp(sprintf(' snap10: %.3f  (%3.2f sec)\n', snap10norm_inf(r), t))
        tic, dmdnorm_inf(r) = norm(chanmin - chan_dmdcts,inf); t=toc;
    disp(sprintf(' dmd: %.3f  (%3.2f sec)\n', dmdnorm_inf(r), t))
    tic, dmdbtnorm_inf(r) = norm(chanmin - chan_btdmd,inf); t=toc;
    disp(sprintf(' dmdbt: %.3f  (%3.2f sec)\n', dmdbtnorm_inf(r), t))
    tic, erabtnorm_inf(r) = norm(chanmin - chan_btera,inf); t=toc;
    disp(sprintf(' erabt: %.3f  (%3.2f sec)\n', erabtnorm_inf(r), t))
    tic, eranorm_inf(r) = norm(chanmin - chan_eracts,inf); t=toc;
    disp(sprintf(' era: %.3f  (%3.2f sec)\n', eranorm_inf(r), t))
    tic, era10norm_inf(r) = norm(chanmin - chan_eracts10,inf); t=toc;
    disp(sprintf(' era10mode: %.3f  (%3.2f sec)\n', era10norm_inf(r), t))
    tic, podmodnorm_inf(r) = norm(chanmin - chan_podmod,inf); t=toc;
    disp(sprintf(' podmod: %.3f  (%3.2f sec)\n', podmodnorm_inf(r), t))
    end
end
%end
%%
ymin = impulse(chanmin,tp);
ymin = ymin';

%% Results for models containing different numbers of modes

markersize = 10;
linewidth = 1.5;
fontsize = 16;
normdata = zeros(length(tp),1);
normmin = zeros(length(tp),1);
normbt = zeros(length(tp),1);
normbpod5 = zeros(length(tp),1);
normbpod10 = zeros(length(tp),1);
normpod = zeros(length(tp),1);
normdmd = zeros(length(tp),1);
normdmdbt = zeros(length(tp),1);
normera = zeros(length(tp),1);
normpodmod = zeros(length(tp),1);

nModes = [3 5 10];
nModes = 3;
for nn = 1:length(nModes);
    
    ybt = impulse(models(nModes(nn)).chan_bt,tp);
    ybt = ybt';
    ybpod10 = impulse(models(nModes(nn)).chan_snap10,tp);
    ybpod10 = ybpod10';
    %ybpod5 = impulse(chan_snap10,tp);
    %ybpod5 = ybpod5';
    ypod = impulse(models(nModes(nn)).chan_pod,tp);
    ypod = ypod';
    
    
    ydmd = impulse(models(nModes(nn)).chan_dmd,tp);
    ydmd = ydmd';
    ydmdbt = impulse(models(nModes(nn)).chan_btdmd,tp);
    ydmdbt = ydmdbt';
    yera = impulse(models(nModes(nn)).chan_era,tp);
    yera = yera';
    ypodmod = impulse(models(nModes(nn)).chan_podmod,tp);
    ypodmod = ypodmod';
    % Plot vorticity fields:
    %{
    tSteps = round(33/dt+1);
    for jj = 1:length(tSteps)
        tStep = tSteps(jj);
        figure
        subplot(3,2,1)
        plotvarinterp(z,y,primal((end/2+1:end),tStep),6,100,1)
        title('Full simulation')
        subplot(3,2,2)
        plotvarinterp(z,y,ybt(end/2+1:end,tStep),6,100,1)
        title('Balanced truncation')
        subplot(3,2,3)
        plotvarinterp(z,y,ypod(end/2+1:end,tStep),6,100,1)
        title('POD projection')
        subplot(3,2,4)
        %plotvarinterp(z,y,ydmd((end/2+1:end),tStep),6,100,1)
        %plotvarinterp(z,y,Bmat(end/2+1:end),6,100,1)
        plotvarinterp(z,y,ybpod10((end/2+1:end),tStep),6,100,1)
        title('Balanced POD projection')
        subplot(3,2,5)
        plotvarinterp(z,y,ydmd((end/2+1:end),tStep),6,100,1)
        title('DMD, POD projection')
        subplot(3,2,6)
        plotvarinterp(z,y,ydmdbt((end/2+1:end),tStep),6,100,1)
        title('DMD, balanced truncation')
        suptitle([num2str(nModes(nn)),' Modes, T = ', num2str(tp(tStep))])
       % if plotRes
    end
    %}
    
%end
%

%for nn = 1:length(nModes);
    for tt = 1:length(tp)
        normdata(tt) = norm(primal(:,tt));
        normmin(tt) = norm(ymin(:,tt));
        normbt(tt) = norm(ybt(:,tt));
        %normbpod5(tt) = norm(ybpod5(:,tt));
        normbpod10(tt) = norm(ybpod10(:,tt));
        normpod(tt) = norm(ypod(:,tt));
        normdmd(tt) = norm(ydmd(:,tt));
        normdmdbt(tt) = norm(ydmdbt(:,tt));
        normera(tt) = norm(yera(:,tt));
        normpodmod(tt) = norm(ypodmod(:,tt));
    end
    %
    
    Colors = get(0,'DefaultAxesColorOrder');
    figure
    %{
    plot(tp,normdata.^2,'LineWidth',linewidth)
    hold on
    %plot(tp,normmin)
    plot(tp,normbt.^2,'LineWidth',linewidth)
    %plot(tp,snap5norm_2)
    %plot(tp,normbpod5)
    plot(tp,normbpod10.^2,'LineWidth',linewidth)
    plot(tp,normpod.^2,'LineWidth',linewidth)
    plot(tp,normdmd.^2,'LineWidth',linewidth)
    plot(tp,normdmdbt.^2,'LineWidth',linewidth)
    plot(tp,normera.^2,'LineWidth',linewidth)
    %}
    plot(tp,normdata.^1,'k','LineWidth',linewidth)
    hold on
    plot(tp,normbt.^1,'--','LineWidth',linewidth,'color',Colors(2,:))
  %  plot(tp,normbpod10.^1,'LineWidth',linewidth)
    plot(tp,normpod.^1,'LineWidth',linewidth,'color',Colors(1,:))
    plot(tp,normdmd.^1,'LineWidth',linewidth,'color',Colors(5,:))
  %  plot(tp,normdmdbt.^1,':','LineWidth',linewidth)
    plot(tp,normera.^1,'LineWidth',linewidth,'color',Colors(6,:))
    
    
     %   plot(tp,normpodmod.^2,'m--','LineWidth',linewidth)

   % plot(tp,normera)
   % legend('Full simulation','Balanced truncation',...
   %     'BPOD model','POD model','DMD model, POD truncation','DMD model, balanced truncation','ERA model, POD truncation')
    legend('Full simulation','BT, BPOD, BTDMD, BTERA models',...
        'POD model','DMD model, POD truncation','ERA model, POD truncation')

    title([num2str(nModes(nn)),' Mode Models'])
    set(gca,'FontSize',fontsize)
    xlabel('Time')
    ylabel('Energy')
    ylabel('$\|(v,\eta)\|_2$','interpreter','latex')
set(gca,'FontSize',fontsize)
    
end

%%
%{
ybt = impulse(chan_bt,tp);
ybt = ybt';
ybpod10 = impulse(chan_snap10,tp);
ybpod10 = ybpod10';
%ybpod5 = impulse(chan_snap5,tp);
%ybpod5 = ybpod5';
ypod = impulse(chan_pod,tp);
ypod = ypod';
ymin = impulse(chan,tp);
ymin = ymin';

ydmd = impulse(chan_dmdcts,tp);
ydmd = ydmd';
ydmdbt = impulse(chan_btdmd,tp);
ydmdbt = ydmdbt';
%}
%%
%{
figure
tStep = 20;
subplot(5,1,1)
plotvarinterp(z,y,primal((end/2+1:end),tStep),6,100,1)
subplot(5,1,2)
%plotvarinterp(z,y,ybt((end/2+1:end),tStep),6,100,1)
%temp = T_pod(:,1:r)*(chan_pod.b);
plotvarinterp(z,y,ymin(end/2+1:end,tStep),6,100,1)
subplot(5,1,3)
plotvarinterp(z,y,ybt(end/2+1:end,tStep),6,100,1)
%temp = chan.c*(chan.b);
%temp = chan.b;
%plotvarinterp(z,y,temp(end/2+1:end),6,100,1)
subplot(5,1,4)
%plotvarinterp(z,y,ydmd((end/2+1:end),tStep),6,100,1)
%plotvarinterp(z,y,Bmat(end/2+1:end),6,100,1)
plotvarinterp(z,y,ypod((end/2+1:end),tStep),6,100,1)
subplot(5,1,5)
plotvarinterp(z,y,ydmdbt((end/2+1:end),tStep),6,100,1)
%}
%% Plot Energies over time


%note that for the optimal IC, hasn't even decayed by 200 seconds, optimal
%at 33 secs (for Re = 1000)
%neither has original initial condition
% For cases where the results are very accurate, the discretization of the
% dmd model contributes the most to the error
%%
% 2-norm of reduced models

figure
ind = 1:rMax;
% h = semilogy(ind, btnorm_2/channorm_2, 'bx-',...
%     ind, snap10norm_2/channorm_2, 'ms--',...
%     ind, podnorm_2/channorm_2, 'r^-',...
%     ind, dmdnorm_2/channorm_2, 'ko-',...
%     ind, dmdbtnorm_2/channorm_2, 'bo--',...  % ind, era10norm_2/channorm_2, 'mx-',... 
%     ind, eranorm_2/channorm_2, 'ro--');%,...   
%     %ind, erabtnorm_2/channorm_2, 'bs--');
h = semilogy(ind, btnorm_2/channorm_2, 'x-',...
    ind, snap10norm_2/channorm_2, 's--',...
    ind, podnorm_2/channorm_2, '^-',...
    ind, dmdnorm_2/channorm_2, 'o-',...
    ind, dmdbtnorm_2/channorm_2, 'x--',...  % ind, era10norm_2/channorm_2, 'mx-',... 
    ind, eranorm_2/channorm_2, 'o--');
set(h,'MarkerSize',markersize, 'LineWidth',linewidth)
set(gca,'FontSize',fontsize)
basename = 'errnorm_2';
leg = legend('Balanced truncation','BPOD model,','POD model',...
    'DMD model, POD truncation','DMD model, balanced truncation',...%'ERA10 modelPOD truncation',...
        'ERA model, POD truncation')%,'ERA model, balanced truncation')
set(leg,'location','southwest')

xlabel('Order of reduced model')
ylabel('$\frac{\|G_r - G\|_2}{\|G\|_2}$','interpreter','latex')
%% inf norms
figure
ind = 1:rMax;
h = semilogy(ind, btnorm_inf/channorm_inf, 'kx-','MarkerSize',markersize, 'LineWidth',linewidth)
hold on
    semilogy(ind, snap10norm_inf/channorm_inf, 's--','color',Colors(2,:),'MarkerSize',markersize, 'LineWidth',linewidth)
    semilogy(ind, podnorm_inf/channorm_inf, '^-','color',Colors(1,:),'MarkerSize',markersize, 'LineWidth',linewidth)
    semilogy(ind, dmdnorm_inf/channorm_inf, 'o-','color',Colors(5,:),'MarkerSize',markersize, 'LineWidth',linewidth)
    semilogy(ind, dmdbtnorm_inf/channorm_inf, '^--','color',Colors(4,:),'MarkerSize',markersize, 'LineWidth',linewidth)  % ind, era10norm_inf/channorm_inf, 'mx-',... 
    semilogy(ind, eranorm_inf/channorm_inf, 's--','color',Colors(6,:),'MarkerSize',markersize, 'LineWidth',linewidth);%,...   
    %ind, erabtnorm_inf/channorm_inf, 'bs--');

set(h,'MarkerSize',markersize, 'LineWidth',linewidth)
set(gca,'FontSize',fontsize)
basename = 'errnorm_2';
leg = legend('Balanced truncation','BPOD model,','POD model',...
    'DMD model, POD truncation','DMD model, balanced truncation',...%'ERA10 modelPOD truncation',...
        'ERA model, POD truncation')%,'ERA model, balanced truncation')
set(leg,'location','northeast')

xlabel('Order of reduced model')
ylabel('$\frac{\|G_r - G\|_{\infty}}{\|G\|_{\infty}}$','interpreter','latex')
%%
figure


plotvarinterp(z,y,Bmat((end/2+1:end)),6,100,1)
xlabel('z')
ylabel('y')
%%
figure
set(gca,'FontSize',fontsize)
plot(cumsum(sig_pod/sum(sig_pod)),'kx-','linewidth',1.5)
xlim([1,10])
xlabel('POD mode number')
ylabel('Cumulative energy fraction')
set(gca,'FontSize',fontsize)
%% plot first few modes for a range of techniques

for modeind = 1:5;
    figure
subplot(3,1,1)
plotvarinterp(z,y,T_pod(end/2+1:end,modeind),6,100,1)
title(['POD mode ', num2str(modeind)])
subplot(3,1,2)
plotvarinterp(z,y,T_snap10(end/2+1:end,modeind),6,100,1)
title(['BPOD mode ', num2str(modeind)])
subplot(3,1,3)
plotvarinterp(z,y,Tinv_snap10(modeind,end/2+1:end),6,100,1)
title(['Adjoint BPOD mode ', num2str(modeind)])
set(gcf,'position',[200 200 400 500])
end
%% plot vel
for modeind = 1:5;
    figure
subplot(3,1,1)
plotvarinterp(z,y,T_pod(1:end/2,modeind),6,100,1)
title(['POD mode ', num2str(modeind)])
subplot(3,1,2)
plotvarinterp(z,y,T_snap10(1:end/2,modeind),6,100,1)
title(['BPOD mode ', num2str(modeind)])
subplot(3,1,3)
plotvarinterp(z,y,Tinv_snap10(modeind,1:end/2),6,100,1)
title(['Adjoint BPOD mode ', num2str(modeind)])
set(gcf,'position',[200 200 400 500])
end
%% test the observability of POD modes

for nMode = 1:10
    ObBPOD(nMode) = T_snap10(:,nMode)'*Wo*T_snap10(:,nMode);
    ObPOD(nMode) = T_pod(:,nMode)'*Wo*T_pod(:,nMode);
    CtPOD(nMode) = T_pod(:,nMode)'*Wc*T_pod(:,nMode);
    Test(nMode) = T_pod(:,nMode)'*Wo^0.5*Wc*Wo^0.5*T_pod(:,nMode);
    
end
figure
semilogy(ObPOD,'linewidth',1.5)
set(gca,'fontsize',14)
%hold on
%plot(ObBPOD)
xlabel('Mode Number')
ylabel('Observability')
%%
figure
semilogy(abs(Test),'linewidth',1.5)
set(gca,'fontsize',14)
%%
figure
semilogy(CtPOD,'linewidth',1.5)
hold on
semilogy(sig_pod(1:10),'linewidth',1.5)
set(gca,'fontsize',14)

xlabel('Mode Number')
ylabel('Controllability')
    
%% Try 3 mode POD models with other modes

modes2useA = [1 3 5 ]
chan_podA    = reduced_model(Amat, Bmat, Cmat, T_pod(:,modes2useA),    T_pod(:,modes2useA)', 3);
modes2useB = [1 2 3 ]
chan_podB    = reduced_model(Amat, Bmat, Cmat, T_pod(:,modes2useB),    T_pod(:,modes2useB)', 3);

%H2 errors
%  A:  pod: 136.667  (0.03 sec)
  %B:   pod: 670.899  (0.02 sec)
 %Hinf errors
 %    pod: 1271.303  (9.53 sec)
 %   pod: 89435.931  (10.30 sec) 
tic, podnorm_A(r) = norm(chanmin - chan_podA,2); t=toc;
    disp(sprintf('    pod: %.3f  (%3.2f sec)', podnorm_A(r), t))
ypodA = impulse(chan_podA,tp);
ypodA = ypodA';
tic, podnorm_B(r) = norm(chanmin - chan_podB,2); t=toc;
 disp(sprintf('    pod: %.3f  (%3.2f sec)', podnorm_B(r), t))
ypodB = impulse(chan_podB,tp);
ypodB = ypodB';
    
    for tt = 1:length(tp)
        normpodA(tt) = norm(ypodA(:,tt));
        normpodB(tt) = norm(ypodB(:,tt));
    end
    
    figure
    plot(tp,normdata,'k','LineWidth',linewidth)
    hold on
    plot(tp,normpodB,'LineWidth',linewidth)
    plot(tp,normpodA,'LineWidth',linewidth)
    legend('True data','POD model, modes 1,2,3','POD model, modes 1,3,5')
    xlabel('Time')
    ylabel('Energy')
set(gca,'FontSize',fontsize)
    %%
    
    
    
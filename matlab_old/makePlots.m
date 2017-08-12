Colors = get(0,'DefaultAxesColorOrder');
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
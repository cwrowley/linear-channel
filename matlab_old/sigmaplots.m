%sigmaplots

load WorkspaceRe1000Tgrowth33dt0p05T1000Ny32.mat
Colors = get(0,'DefaultAxesColorOrder');
%%
figure
 sigmaplot(models(10).chan_bt)
 
 %%
nModes = [3 5 10];
nModes = 3;

for nn = 1:length(nModes);
    figure
     p1 = sigmaplot(models(nModes(nn)).chan_bt,'b')
     L = findobj(gcf,'type','line')
     set(L,'Color' ,[ 0 0 0]);
    % set(p1,'cor',[ 0 0 0])
     hold on
     sigmaplot(models(nModes(nn)).chan_snap10,'b')%,'--')
      L2 = findobj(gcf,'type','line','-and','Color','b');
     set(L2,'Color' ,Colors(2,:));
     sigmaplot(models(nModes(nn)).chan_pod,'b')
     L3 = findobj(gcf,'type','line','-and','Color','b');
     set(L3,'Color' ,Colors(1,:));
     sigmaplot(models(nModes(nn)).chan_dmd,'b')
     L4 = findobj(gcf,'type','line','-and','Color','b');
     set(L4,'Color' ,Colors(5,:));
     sigmaplot(models(nModes(nn)).chan_btdmd,'b:')%,':')
     L5 = findobj(gcf,'type','line','-and','Color','b');
     set(L5,'Color' ,Colors(4,:));
     sigmaplot(models(nModes(nn)).chan_era,'b')
     L6 = findobj(gcf,'type','line','-and','Color','b');
     set(L6,'Color' ,Colors(6,:));
     set(gca,'fontsize',14)
     grid on
     xlim([1e-4,10])
    leg = legend('Balanced truncation',...
        'BPOD model','POD model','DMD model, POD truncation','DMD model, balanced truncation','ERA model, POD truncation')
    set(leg,'location','southwest')
    title([num2str(nModes(nn)),' Mode Models'],'fontsize',16)
    %print(
end

%%
nModes = 3;Lwidth = 1.5;
%w = logspace(-4,1,100);
for nn = 1:length(nModes);
    [sv,w] = sigma(models(nModes(nn)).chan_bt);
    figure
    loglog(w,sv,'k')
    hold on
    [sv,w] = sigma(models(nModes(nn)).chan_snap10);
    loglog(w,sv,'--','color',Colors(2,:),'linewidth',Lwidth)
    [sv,w] = sigma(models(nModes(nn)).chan_pod);
    loglog(w,sv,'Color' ,Colors(1,:),'linewidth',Lwidth)
    [sv,w] = sigma(models(nModes(nn)).chan_dmd);
    loglog(w,sv,'Color' ,Colors(5,:),'linewidth',Lwidth)
    [sv,w] = sigma(models(nModes(nn)).chan_btdmd);
    loglog(w,sv,':','Color' ,Colors(4,:),'linewidth',Lwidth)
    [sv,w] = sigma(models(nModes(nn)).chan_era);
    loglog(w,sv,'Color' ,Colors(6,:),'linewidth',Lwidth)
    
    set(gca,'fontsize',14)
     grid on
     xlim([1e-3,1])
    xlabel('Frequency (rad/s)')
    ylabel('Singular Values (dB)')
    leg = legend('Balanced truncation',...
        'BPOD model','POD model','DMD model, POD truncation','DMD model, balanced truncation','ERA model, POD truncation')
    set(leg,'location','southwest')
    title([num2str(nModes(nn)),' Mode Models'],'fontsize',16)
end
%%
  nModes = 3;Lwidth = 1.5;
%w = logspace(-4,1,100);
for nn = 1:length(nModes);
    %[sv,w] = sigma(models(nModes(nn)).chan_bt);
    figure
    %loglog(w,sv,'k')
  %  hold on
    [sv,w] = sigma(models(nModes(nn)).chan_snap10);
    loglog(w,sv,'-','color',Colors(2,:),'linewidth',Lwidth)
    hold on
    [sv,w] = sigma(models(nModes(nn)).chan_pod);
    loglog(w,sv,'Color' ,Colors(1,:),'linewidth',Lwidth)
    [sv,w] = sigma(models(nModes(nn)).chan_dmd);
    loglog(w,sv,'Color' ,Colors(5,:),'linewidth',Lwidth)
    %[sv,w] = sigma(models(nModes(nn)).chan_btdmd);
    %loglog(w,sv,':','Color' ,Colors(4,:),'linewidth',Lwidth)
    [sv,w] = sigma(models(nModes(nn)).chan_era);
    loglog(w,sv,'Color' ,Colors(6,:),'linewidth',Lwidth)
    
    set(gca,'fontsize',14)
     grid on
     xlim([1e-3,1])
    xlabel('Frequency (rad/s)')
    ylabel('Singular Values (dB)')
    leg = legend('BT, BPOD, BTDMD, BTERA models',...
        'POD model','DMD model, POD truncation','ERA model, POD truncation')
    set(leg,'location','southwest')
    title([num2str(nModes(nn)),' Mode Models'],'fontsize',16)
end  

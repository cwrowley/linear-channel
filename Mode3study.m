nModes = [3];

for nn = 1:length(nModes);
    
    [ypod,~,xpod] = impulse(models(nModes(nn)).chan_pod,tp);
    ypod = ypod';


%for nn = 1:length(nModes);
    for tt = 1:length(tp)

        normpod(tt) = norm(ypod(:,tt));
        normdmd(tt) = norm(ydmd(:,tt));

    end
    %
    figure
    plot(tp,normdata,'LineWidth',linewidth)
    hold on

    plot(tp,normpod,'LineWidth',linewidth)
    legend('Full simulation','POD model')
    title([num2str(nModes(nn)),' Mode Models'])
    set(gca,'FontSize',fontsize)
    xlabel('Time')
    ylabel('Energy')
set(gca,'FontSize',fontsize)
    
end
%%
figure
plot(tp,xpod');
hold on
plot(tp,normpod,'LineWidth',linewidth)
legend('x_1','x_2','x_3','Energy')
set(gca,'FontSize',fontsize)
    xlabel('Time')
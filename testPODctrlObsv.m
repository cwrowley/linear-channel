% Analyze controllability and observabitily

load WorkspaceRe1000Tgrowth33dt0p05.mat
%%
for nMode = 1:10
    ObBPOD(nMode) = T_snap10(:,nMode)'*Wo*T_snap10(:,nMode);
    ObPOD(nMode) = T_pod(:,nMode)'*Wo*T_pod(:,nMode);
    CtPOD(nMode) = T_pod(:,nMode)'*Wc*T_pod(:,nMode);
    Test(nMode) = T_pod(:,nMode)'*Wc^0.5*Wo*Wc^0.5*T_pod(:,nMode);
    Test2(nMode) = T_pod(:,nMode)'*Wo^0.5*Wc*Wo^0.5*T_pod(:,nMode);
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
hold on
%semilogy(abs(Test2),'linewidth',1.5)
legend('|u_i^*Wc^{1/2}WoWc^{1/2}u_i','u_i^*Wo^{1/2}WcWo^{1/2}u_i|')
set(gca,'fontsize',14)
xlabel('Mode Number')
%%
figure
semilogy(CtPOD/dt,'linewidth',1.5)
hold on
semilogy(sig_pod(1:10).^2,'linewidth',1.5)
set(gca,'fontsize',14)
legend('Mode Controllability','POD Eigenvalue')
xlabel('Mode Number')
%ylabel('Controllability')
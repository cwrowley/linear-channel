% check subspaces

%load WorkspaceRe1000Tgrowth33dt0p05.mat
load subspaces.mat
%% Compute angles between subspaces
rMax = 8;
for r = 1:rMax
    theta_pod_bt(r) = subspace(T_pod(:,1:r),T(:,1:r));
    thetapod_bpod(r) = subspace(T_pod(:,1:r),T_snap10(:,1:r));
    theta_pod_adjbt(r) = subspace(T_pod(:,1:r),Tinv(1:r,:)');
    theta_pod_adjbpod(r) = subspace(T_pod(:,1:r),Tinv_snap10(1:r,:)');
    
    theta_bt_bpod(r) = subspace(T(:,1:r),T_snap10(:,1:r));
    theta_bt_btadj(r) = subspace(T(:,1:r),Tinv(:,1:r));
end


figure
plot(1:rMax,theta_pod_bt)
hold on
plot(1:rMax,thetapod_bpod)
plot(1:rMax,theta_pod_adjbt)
plot(1:rMax,theta_pod_adjbpod,'--')
plot(1:rMax,theta_bt_bpod,'-')
plot(1:rMax,theta_bt_btadj,'--')
set(gca,'FontSize',14)
legend('POD,BT','POD,BPOD','POD,BTadj','POD,BPODadj','BT,BPOD','BT,BTadj')
xlabel('Subspace dimension')
ylabel('\theta')

%% Look at Grassmanian metric

for r = 1:rMax
    d_pod_bt(r) = subspace_metric(T_pod(:,1:r),T(:,1:r));
    d_pod_bpod(r) = subspace_metric(T_pod(:,1:r),T_snap10(:,1:r));
    d_pod_adjbt(r) = subspace_metric(T_pod(:,1:r),Tinv(1:r,:)');
    d_pod_adjbpod(r) = subspace_metric(T_pod(:,1:r),Tinv_snap10(1:r,:)');
    d_bt_bpod(r) = subspace_metric(T(:,1:r),T_snap10(:,1:r));
    d_bt_adjbt(r) = subspace_metric(T(:,1:r),Tinv(:,1:r));
end
figure
plot(1:rMax,d_pod_bt)
hold on
plot(1:rMax,d_pod_bpod,'--')
plot(1:rMax,d_pod_adjbt)
plot(1:rMax,d_pod_adjbpod,'--')
plot(1:rMax,d_bt_bpod,'--')
plot(1:rMax,d_bt_adjbt,'--')
set(gca,'FontSize',14)
legend('POD,BT','POD,BPOD','POD,BTadj','POD,BPODadj','BT,BPOD','BT,BTadj')
xlabel('Subspace dimension')
ylabel('Distance')
%resolution studies
clear all
close all
Nzs = [8 16 32 ];
Nys = [8 16 32 ];
for NyInd = 1:length(Nzs);
    for NzInd = 1:length(Nys);
        Nz = Nzs(NzInd);
        Ny = Nys(NyInd);
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
        
        
        %Aadj is now defined in consistant coords
        [Amat, Aadj, Bmat, Mmat, Npts, y, z] = define_eqns(Ny, Nz, Re);
        Cmat = eye(2*Npts);
        %%
        if icFlag
            Tgrowths =[30:0.1:36];%[1 30 50 200];%[0.1:0.1:100]; %optimal time
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
        %  disp('Computing minimal realization')
        %  tic
        %  chanmin = minreal(chan); % use for computing norms
        % t = toc; disp(sprintf('    time = %.2f',t));
        %    end
        
        [T_pod, sig_pod, primal] = compute_pod(chan,tp);
        
        ResolutionTest(NyInd,NzInd).primal = primal;
        ResolutionTest(NyInd,NzInd).T_pod = T_pod;
        ResolutionTest(NyInd,NzInd).sig_pod = sig_pod;
        ResolutionTest(NyInd,NzInd).y = y;
        ResolutionTest(NyInd,NzInd).z = z;
        ResolutionTest(NyInd,NzInd).growthInd = growthInd;
        %  ResolutionTest(NyInd).chanmin = chanmin;
    end
end
%%

figure
for NyInd = 1:3;
    semilogy(     ResolutionTest(NyInd).sig_pod)
    hold on
end

%%
%norm(ResolutionTest(3).chanmin - ResolutionTest(1).chanmin,2);
%norm(ResolutionTest(3).chanmin - ResolutionTest(2).chanmin,2);

%%
%ymin = impulse(chanmin,tp);
%ymin = ymin';
for tt = 1:length(tp)
normdataT(tt) = norm(ResolutionTest(3,3).primal(:,tt));
end
Error = zeros(3,3);
TotalE = zeros(3,3);
%% All
figure
for NzInd = 1:3;
    for NyInd = 1:3;
        for tt = 1:length(tp)
           
            normdata(tt) = norm(ResolutionTest(NyInd,NzInd).primal(:,tt));
            %normmin(tt) = norm(ymin(:,tt));
        end

          plot(tp,normdata)
          hold on
          
 Error(NzInd,NyInd) = norm(normdata-normdataT);
 
 TotalE(NzInd,NyInd) = sum(normdata)*dt;
    end
end
    
legend('Nz 8, Ny 8','Nz 8, Ny 16','Nz 8, Ny 32','Nz 16, Ny 8','Nz 16, Ny 16',...
    'Nz 16, Ny 32','Nz 32, Ny 8','Nz 32, Ny 16','Nz 32, Ny 32')

%%
%% Change Ny
fontsize = 16;
figure
for NzInd = 2;
    for NyInd = 1:3;
        for tt = 1:length(tp)
           
            normdata(tt) = norm(ResolutionTest(NyInd,NzInd).primal(:,tt));
            %normmin(tt) = norm(ymin(:,tt));
        end

          plot(tp,normdata.^2)
          hold on
          
 Error(NzInd,NyInd) = norm(normdata-normdataT);
 
 TotalE(NzInd,NyInd) = sum(normdata)*dt;
    end
end
    
legend('Nz 16, Ny 8','Nz 16, Ny 16',...
    'Nz 16, Ny 32','Nz 16, Ny 64')
    xlabel('Time')
    ylabel('Energy')
set(gca,'FontSize',fontsize)
title('Changing Ny')
%% Change Nz
figure
for NzInd = 1:3;
    for NyInd = 2;
        for tt = 1:length(tp)
           
            normdata(tt) = norm(ResolutionTest(NyInd,NzInd).primal(:,tt));
            %normmin(tt) = norm(ymin(:,tt));
        end

          plot(tp,normdata.^2)
          hold on
          
 Error(NzInd,NyInd) = norm(normdata-normdataT);
 
 TotalE(NzInd,NyInd) = sum(normdata)*dt;
    end
end
    
legend('Nz 8, Ny 16','Nz 16, Ny 16',...
    'Nz 32, Ny 16')
    xlabel('Time')
    ylabel('Energy')
set(gca,'FontSize',fontsize)
title('Changing Nz')

%%
    figure
    tind = 1;
subplot(3,1,1)
plotvarinterp(ResolutionTest(2,2).z,ResolutionTest(2,2).y,ResolutionTest(2,2).primal(end/2+1:end,tind),6,100,1)
%title(['POD mode ', num2str(modeind)])
subplot(3,1,2)
plotvarinterp(ResolutionTest(3,3).z,ResolutionTest(3,3).y,ResolutionTest(3,3).primal(end/2+1:end,tind),6,100,1)
%title(['POD mode ', num2str(modeind)])
subplot(3,1,3)
plotvarinterp(ResolutionTest(3,2).z,ResolutionTest(3,2).y,ResolutionTest(3,2).primal(end/2+1:end,tind),6,100,1)
%title(['POD mode ', num2str(modeind)])


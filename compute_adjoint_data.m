function [dualdat] = compute_adjoint_data(Amat, Phi, tp)

%Amat is adjoint matrix!

N = size(Amat,1);
projrank = size(Phi,2);
nt = length(tp);



%dualsys = ss(Amat', Phi, eye(N), 0);
dualsys = ss(Amat, Phi, eye(N), 0);


% run simulation for each POD mode
%  POD modes contained in columns of Phi
dualdat = zeros(nt*projrank, N);
for nsim=1:projrank
    disp(sprintf('  Integrating dual problem %d of %d...',nsim,projrank))
    x = initial(dualsys,Phi(:,nsim),tp); %Phi is now in vel,vort coords
    offset = (nsim-1)*nt;
    dualdat(offset+1:offset+nt,:) = x;
end
dualdat = dualdat';
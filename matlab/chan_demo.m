% Script demonstrating object-oriented version of linear channel code
% Clancy Rowley, 12 Aug 2017

ny = 32;
nz = 16;
Re = 1000;

linchan = LinearChannel(ny, nz, Re);

t_opt = 32.9;
B = linchan.optimal_perturbation(t_opt);

[vel,vort] = linchan.vec2field(B);

% Plot optimal perturbation (initial condition)
figure(1)
subplot(211)
linchan.plotfield(vel)
colorbar
title('Wall-normal velocity')

subplot(212)
linchan.plotfield(vort)
colorbar
title('Wall-normal vorticity')

% Construct a linear system and simulate
C = eye(2*linchan.npts);
sys = ss(linchan.A,B,C,0);
t_final = 300;
[y,t] = impulse(sys, t_final);

% Plot 'energy' (not the true kinetic energy)
energy = dot(y', y');
figure(2)
plot(t, energy);
xlabel('Time')
ylabel('Energy |x|^2')

% Plot solution at a later time
ind = find(t > t_opt, 1);
[vel, vort] = linchan.vec2field(y(ind,:));
figure(3)
subplot(211)
linchan.plotfield(vel)
colorbar
title(sprintf('Wall-normal velocity, t = %f', t(ind)))
subplot(212)
linchan.plotfield(vort)
colorbar
title(sprintf('Wall-normal vorticity, t = %f', t(ind)))

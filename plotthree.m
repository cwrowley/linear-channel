function plotthree(z,y,one,two,three,ind,signum);
% v1 eta1
% v2 eta2
% v3 eta3

clev = 14;
Npts = size(one,2)/2;
lrmarg = 0.05;
tbmarg = 0.05;
midv = 0.05;
midh = 0;

wid = (1-2*lrmarg-midv)/2;
ht = (1-2*tbmarg-2*midh)/3;

subplot('position',[lrmarg 1-tbmarg-ht wid ht])
plotvar(z,y,signum(1)*one(ind,1:Npts),clev);
axis([0 z(end) y(end) y(1)]);
%axis([0 2 y(end) y(1)]);
%title('Normal velocity')

%subplot(322)
subplot('position',[0.5+midv/2 1-tbmarg-ht wid ht])
plotvar(z,y,signum(1)*one(ind,Npts+1:end),clev);
axis([0 z(end) y(end) y(1)]);
%axis([0 2 y(end) y(1)]);
%title('Normal vorticity')

%subplot(323)
subplot('position',[lrmarg tbmarg+ht+midh wid ht])
plotvar(z,y,signum(2)*two(ind,1:Npts),clev);
axis([0 z(end) y(end) y(1)]);
%axis([0 2 y(end) y(1)]);
%subplot(324)
subplot('position',[0.5+midv/2 tbmarg+ht+midh wid ht])
plotvar(z,y,signum(2)*two(ind,Npts+1:end),clev);
axis([0 z(end) y(end) y(1)]);
%axis([0 2 y(end) y(1)]);

%subplot(325)
subplot('position',[lrmarg tbmarg wid ht])
plotvar(z,y,signum(3)*three(ind,1:Npts),clev);
axis([0 z(end) y(end) y(1)]);
%axis([0 2 y(end) y(1)]);
%subplot(326)
subplot('position',[0.5+midv/2 tbmarg wid ht])
plotvar(z,y,signum(3)*three(ind,Npts+1:end),clev);
axis([0 z(end) y(end) y(1)]);
%axis([0 2 y(end) y(1)]);

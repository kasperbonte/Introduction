clearvars
clc
close all
%%
% Analytical and numerical solution wild boars
N(1) = 5*10^3;
r = 0.48;
dt = 1;
totTime = 10;
time = dt:dt:totTime;
tic
for t = 2:numel(time)
        dN = ((r*N(t-1))*dt);
        N(t) = N(t-1)+dN;
end
toc
Popan(1) = N(1);
Popan(1,2:numel(time)) = N(1)*exp(r*time(1,1:numel(time)-1));
figure
plot(time,N)
hold on
plot(time,Popan)
%%
% Carrying capicity wild boars numerically
clearvars
clc
close all
area = 30528
density = 20
K = density * area
r = 0.48 
N(1) = 5000
dt = 1
totTime = 100
time = dt:dt:totTime
for t = 2:numel(time)
    dN = r*((K-N(t-1))/K)*N(t-1)
    N(t) = N(t-1)+dN
end
plot(time,N)
%%
% analytical solution carrying capacity 
clearvars
clc
close all
area = 30528
density = 20
K = density * area
r = 0.48 
N(1) = 5000
dt = 1
totTime = 100
time = dt:dt:totTime
Popan(1) = N(1)
Popan(1,2:numel(time)) = ((-K*N(1)*exp(r*time(1,1:numel(time)-1)))/(N(1)-K))/(1-(N(1)*exp(r*time(1,1:numel(time)-1)))/(N(1)-K))
%%
%Hunting
clearvars
clc
close all
N0=5000
r=0.48
hunting_effort=0.48+0.48*0.2 %%Stable population if this equals r
dt=1
N=[N0]
n=N0
totTime = 100
time = dt:dt:totTime
H = hunting_effort*N
Popan(1) = N0
for t = 2:numel(time)
  Popan(1,t) = Popan(1,t-1)+(r*Popan(1,t-1)-H)
  H = hunting_effort*Popan(1,t)
end
plot(time,Popan)

%%
%wolves
clearvars
clc
close all
r = 0.48
c = 0.01
d = 0.24
e = 0.005
Nboars = 100
Nwolves = 25
tijd = 1:100
for t = 1:99
    N = Nboars(t)*(1+r)-c*Nwolves(t)*Nboars(t);
    W = Nwolves(t)*(1-d)+e*Nwolves(t)*N
    Nboars = [Nboars,N]
    Nwolves = [Nwolves, W]
    
end
plot(tijd, Nboars)
analytical_error=e*Nboars-d*log(Nboars)-r*log(Nwolves)+c*Nwolves
%%
%Matlab solver
clearvars
clc
close all
r = 0.48;
c = 0.01;
d = 0.24;
e = 0.005;
W = 25
y0 = [0 100]
timeperiod = [0 100]
[t,y] = ode45(@(t,y)odefnc(t,y,r,c,W), timeperiod, initial)

% https://www.youtube.com/watch?v=5r1zFz5CB2w link info method + https://www.youtube.com/watch?v=nTRNrvJmixw
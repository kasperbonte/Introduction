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
tic
r = 0.48
c = 0.01
d = 0.24
e = 0.005
Nboars = 100
Nwolves = 25
dt = 1
Tottime = 100
timeperiod= dt:dt:Tottime
for t = 2:numel(timeperiod)
    N = Nboars(t-1)+((r*Nboars(t-1)-c*Nwolves(t-1)*Nboars(t-1))*dt);
    W = Nwolves(t-1)+((-d*Nwolves(t-1)+e*Nwolves(t-1)*N)*dt)
    Nboars = [Nboars,N]
    Nwolves = [Nwolves, W]
    
end
plot(timeperiod, Nboars,timeperiod,Nwolves)
toc
elapsedTime1 = toc

analytical_error=e*Nboars-d*log(Nboars)-r*log(Nwolves)+c*Nwolves
%%
%Matlab solver
clearvars
clc
close all
tic
r = 0.48;
c = 0.01;
d = 0.24;
e = 0.005;
W = 25;
Initial = [100,25];     % For a differential equation the initial conditions are necessary to solve the equation 
timeperiod = [0 100]; % this is the interval of integration of the differential equation
Functie = @(t,y)[r*y(1)-c*y(2)*y(1);-d*y(2)+e*y(2)*y(1)]; % this commant defines the equations that need to be solved
ode45(Functie,timeperiod,Initial) % specific syntax to solve the equation and plot it
toc
elapsedTime2 = toc
Timediff = elapsedTime2-elapsedTime1

% https://www.youtube.com/watch?v=5r1zFz5CB2w link info method + https://www.youtube.com/watch?v=nTRNrvJmixw

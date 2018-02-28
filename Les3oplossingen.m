clearvars
clc
close all


r=0.48;
c=0.01;
d=0.24;
e=0.005;

Nini = 100;
Wini=25;
 
totT = 15;
N(1) = Nini;
W(1)=Wini;
dt = linspace(0.001,1,100);

% dt=0.01; % gebruik verschillende dt's van 1 => 0.1 en zie het effect op de grafiek => numerische stabiliteit
for p = 1:numel(dt);
    dt = dt(p);
    time=dt:dt:totT;
    N = zeros(1,numel(time));
    W = zeros(1,numel(time));
    for t=2:numel(time);
            
            N(1) = Nini;
            W(1) = Wini;
            dN = (r*N(t-1)-c*W(t-1)*N(t-1))*dt;
            N(t) = N(t-1)+dN;
        
            dW = (-d*W(t-1)+e*W(t-1)*N(t-1))*dt;
            W(t) = W(t-1)+dW;
      
    end
    H_model = Analerror(N,W,d,e,r,c);
    H_true = Analerror(N(1),W(1),d,e,r,c);
    H_RMSE(1,p) = RootMSE(H_true,H_model);
    H_RMSE_P(1,p) = H_RMSE(1,p)/H_true ;
    dt = linspace(0.001,1,100);
        
end

%Figure
% figure
% subplot(2,1,1)
% plot(time,N)
% hold on
% plot(time, W)
% legend('N','W')
% hold off
% xlabel('Time'),ylabel('Population')
% 
% subplot(2,1,2)
% plot(W,N)
% hold off

figure
plot(dt,H_RMSE)
xlabel('timestep')
ylabel('RMSE error')

%% Solve ODE with matlab solver and adapt so the accuracy is similar as own solved numerical equation:
clearvars
clc
close all
r=0.48;
c=0.01;
d=0.24;
e=0.005;
Timeperiod = [0 100]
Initial = [100, 25]
dt = 0.05
Timeperiod =[0:dt:100]
% Function = @(t,y)[y(1)*r-c*y(2)*y(1);y(2)*(-d)+e*y(2)*y(1)]; % https://nl.mathworks.com/help/matlab/matlab_prog/matlab-operators-and-special-characters.html search @ description @
% [t,y] = ode45(Function,Timeperiod,Initial)
const.r = r
const.c = c
const.d = d
const.e = e
Function = @(t,y)calcrate(t,y,const)
options = odeset('RelTol',1e-5) %The tolerances are used to limit the local discretization error: If the difference between a high-accuracy and low-accuracy integration is higher than one of the tolerances, the step size is reduced.
[t,p] = ode45(Function,Timeperiod,Initial,options)
% [t,p] = ode45(@(t,y)calcrate(t,y,const),Timeperiod,Initial)
figure('color','white')
plot(p(:,1),(p(:,2)))
xlabel('Boars')
ylabel('Wolves')

%% Insert human hunting in the system
clearvars
clc
close all

r=0.48;
c=0.01;
d=0.24;
e=0.005;
alpha = 0.1

Nini = 100;
Wini=25;


totT = 200;
dt=0.001;
time=dt:dt:totT;

N(1) = Nini;
W(1)=Wini;
H(1) = 0;

for t=2:numel(time)
%         dH = alpha*N(t-1);
%         H(t) = dH;
%         dN = ((r*N(t-1)-c*W(t-1)*N(t-1))*dt-0.5*dt)-H(t) ´%Insert human hunting with sniper
        dN = (r*N(t-1)-c*W(t-1)*N(t-1))*dt-0.5*dt;% Bij hunting = 0.5 trager naar een diepte punt de boars, de wolves blijven op een constante populatie dan
%         dN = (r*N(t-1)-c*W(t-1)*N(t-1))*dt-1*dt; % Bij hunting = 1 snel
%         naar een diepte punt de boars, de wolves blijven op een constante
%         populatie dan
        N(t) = N(t-1)+dN;
        
        dW = (-d*W(t-1)+e*W(t-1)*N(t-1))*dt;
        W(t) = W(t-1)+dW;
end
%Figure
figure
subplot(2,1,1)
ylim([0 300])
plot(time,N)
hold on
plot(time, W)
ylim([-100 300])
legend('N','W')
hold off
xlabel('Time'),ylabel('Population')

subplot(2,1,2)
plot(W,N)
ylim([-100 300])
xlim([-100,300])

%%
% Insert roe deer + hunting (0.1 per year)
clearvars
clc
close all

r=0.48;
r2 =0.48;
c2 = 0.01;
c=0.01;
d=0.24;
e=0.005;
e2 = 0.005;
Initial = [100,25,80]
dt = 0.01
Timeperiod = [0:dt:100]
Functie = @(t,y)[(r*y(1)-c*y(2)*y(1))-0.1*dt;(-d*y(2)+e*y(2)*y(1)+e2*y(2)*y(3))-0.1*dt;(y(3)*r2-y(3)*y(2)*c2)-0.1*dt]; 
[t,y] = ode45(Functie,Timeperiod,Initial) 
Boars = y(:,1)
Wolves = y(:,2)
Deer = y(:,3)
% plot(t(:,1),Boars,t(:,1),Wolves,t(:,1),Deer)
% legend('Boars','Wolves','Deer')
figure
scatter3(Wolves,Boars,Deer,40,Timeperiod,'filled')
ax = gca
xlabel('Amount of wolves')
ylabel('Amount of boars')
zlabel('Amount of deers')
cb = colorbar;       % https://nl.mathworks.com/help/matlab/examples/visualizing-four-dimensional-data.html
cb.Label.String = 'Timeperiod';















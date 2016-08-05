%%===== SCRIPT FOR TRACKING NETWORK ===================================%%
clear all
clc
close all

c=[-40,-20,0,20,40];
x=[-2,-1,0,1,2];

lambdas=-0; %A si 1D
lambdad=0.005;%0.1;
mu=1; %L2 penality
nu=0.1; %L1 penality

N=4;
Gamma=1*[-2 -1 1 2];
FastConnections=Gamma'*Gamma+mu*lambdad^2*eye(N); %3x3 matrix
SlowConnections=Gamma'*(lambdas+lambdad)*Gamma; %3x3 matrix

V(1)=0;
v(:,1)=zeros(N,1);
x_appr(1)=0;



stimulus=[x(1)*ones(1,24),x(2)*ones(1,25),x(3)*ones(1,25),x(4)*ones(1,25),x(5)*ones(1,25)];
%stimulus=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x(2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x(4),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x(5),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; %on commence avc stimulus constant
myx(1)=x(1); %myx <-- x du papier
t=1;
E(t)=4
r(1,:)=zeros(1,N);
T=1/2*(nu*lambdad+nu*lambdad^2+Gamma.^2);
o(:,1)=zeros(N,1);

while t<=124 %E(t)>0.1
V=Gamma'*(myx(t)-x_appr(t))-mu*lambdad*r(t,:)';
v(:,t)=Gamma'*(myx(t)-x_appr(t))-mu*lambdad*r(t,:)';
o(:,t)=zeros(N,1);
for i=1:N
if V(i)>T(i)
    disp(['At time t=',num2str(t),' neuron ',num2str(i),' spikes.']);
    disp(['V=',num2str(V(i)),' > T=',num2str(T(i)),'.']);
o(i,t)=o(i,t)+1;
end
end
x_appr(t+1)=x_appr(t)-lambdad*x_appr(t)+Gamma*o(:,t);
r(t+1,:)=r(t,:)-lambdad*r(t,:)+lambdad*o(:,t)';
myx(t+1)=myx(t)+lambdas*myx(t)+stimulus(t);
E(t+1)=(myx(t+1)-x_appr(t+1))^2+nu*sum(abs(r(t+1,:)))+mu*sum(r(t+1,:).^2);

t=t+1;
end

figure; 
plot(stimulus); hold on; plot(myx,'g'); plot(x_appr,'r'); hold off
title(['N=',num2str(N),'        Gamma=[',num2str(Gamma),']        lambdas=',num2str(lambdas),'        lambdad=',num2str(lambdad),'        mu=',num2str(mu),'        nu=',num2str(nu)]);

figure;
for i=1:N
subplot(N,1,i)
plot(o(i,:));
title(['Spike train for neuron ',num2str(i),'.']);
end

figure;
imagesc(o);
title('spike trains')
figure;
imagesc(exp(v-min(min(v))));
title('voltage')

figure; imagesc(FastConnections); colorbar; title('Fast connections')
figure; imagesc(SlowConnections); colorbar; title('Slow connections')

figure;
hold on;
plot(v(1,:),'c')
plot(1.8*o(1,:),'b')
plot(T(1)*ones(size(v(1,:))),'r')
hold off


%%===== SCRIPT FOR TRACKING NETWORK in OPTIC TECTUM =====================%%
% clear all
clc
close all

J=4; % dimension du stimulus et de la variable descriptive x

lambdas=-0.9*eye(J); %A si 1D
lambdad=0.61;%0.1;
mu=10; %L2 penality
nu=1; %L1 penality

N=312;
Gamma=zeros(J,N);

nthatpref=zeros(J,1);
for d=1:4
nthatpref(d)=length(find(PrefDir==d));
end

compteur=zeros(4,1);
for n=1:N

    compteur(PrefDir(n))=compteur(PrefDir(n))+1;
    
if compteur(PrefDir(n))<=floor(nthatpref(PrefDir(n))/2)
Gamma(PrefDir(n),n)=10*0.1;
else
Gamma(PrefDir(n),n)=-10*0.1;
end
end

FastConnections=Gamma'*Gamma+mu*lambdad^2*eye(N); %NxN matrix
SlowConnections=Gamma'*(lambdas+lambdad*eye(J))*Gamma; %NxN matrix

V(1)=0;
v(:,1)=zeros(N,1);
x_appr(:,1)=zeros(J,1);

%idee : stimulus taille 4xT binaire, au temps t<T 1 à l'index de la direction présentée, 0 sinon.
%interet : encoder notion de direction preferee d'apres resultats de pop vector decoding.

stimulus=zeros(J,16);
for s=1:J
stimulus(s,(s-1)*J+1:s*J)=ones(1,J);
end
%%

myx(:,1)=stimulus(:,1); %myx <-- x du papier
t=1;
E(t)=0;
r(1,:)=zeros(1,N);
T=1/2*((nu*lambdad+nu*lambdad^2)*ones(N,1)+diag(Gamma'*Gamma));
o(:,1)=zeros(N,1);

while t<=size(stimulus,2) %E(t)>0.1
V=Gamma'*(myx(:,t)-x_appr(:,t))-mu*lambdad*r(t,:)';
v(:,t)=Gamma'*(myx(:,t)-x_appr(:,t))-mu*lambdad*r(t,:)';
o(:,t)=zeros(N,1);
for i=1:N
if V(i)>T(i)
%     disp(['At time t=',num2str(t),' neuron ',num2str(i),' spikes.']);
%     disp(['V=',num2str(V(i)),' > T=',num2str(T(i)),'.']);
o(i,t)=o(i,t)+1;
end
end
x_appr(:,t+1)=x_appr(:,t)-lambdad*x_appr(:,t)+Gamma*o(:,t);
r(t+1,:)=r(t,:)-lambdad*r(t,:)+lambdad*o(:,t)';
myx(:,t+1)=myx(:,t)+lambdas*myx(:,t)+stimulus(:,t);
[~,idx]=max(myx(:,t+1));
[~,idx_appr]=max(x_appr(:,t+1));
for l=1:J
if l==idx
myx(l,t+1)=1;
else
myx(l,t+1)=0;
end
if l==idx_appr
x_appr(l,t+1)=1;
else
x_appr(l,t+1)=0;
end

end
E(t+1)=norm(myx(:,t+1)-x_appr(:,t+1),2)^2+nu*sum(abs(r(t+1,:)))+mu*sum(r(t+1,:).^2);

t=t+1;
end

%% FIG

figure;
plot(E);
title('Error')

figure;
subplot(3,1,1)
imagesc(stimulus)
subplot(3,1,2)
imagesc(myx)
subplot(3,1,3)
imagesc(x_appr)

figure; imagesc(FastConnections); colorbar; title('Fast connections')
figure; imagesc(SlowConnections); colorbar; title('Slow connections')

%%
% 
% % StimulusForPlot=
% figure; 
% plot(stimulus); hold on; plot(myx,'g'); plot(x_appr,'r'); hold off
% title(['N=',num2str(N),'        Gamma=[',num2str(Gamma),']        lambdas=',num2str(lambdas),'        lambdad=',num2str(lambdad),'        mu=',num2str(mu),'        nu=',num2str(nu)]);
% 
% figure;
% for i=1:N
% subplot(N,1,i)
% plot(o(i,:));
% title(['Spike train for neuron ',num2str(i),'.']);
% end
% 
% figure;
% imagesc(o);
% title('spike trains')
% figure;
% imagesc(exp(v-min(min(v))));
% title('voltage')
% 
% figure; imagesc(FastConnections); colorbar; title('Fast connections')
% figure; imagesc(SlowConnections); colorbar; title('Slow connections')
% 
% figure;
% hold on;
% plot(v(1,:),'c');
% plot(1.8*o(1,:),'b');
% plot(T(1)*ones(size(v(1,:))),'r');
% hold off
% 
% figure; 
% plot(E);
% title('error')
% 
% figure;
% hold on
% plot(stimulus,'b');
% plot(x_appr,'c');
% hold off

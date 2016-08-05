%%===========================decoding optic tectum ds====================%%

clear all; 
close all; 
clc
mypath='/users/zfne/mehyaoui/Desktop/phD_Lab2/Donnees Optic Tectum/';
load([mypath,'141024_GCaMP5G_7dpf_LongStim_F1/Results/141024_GCaMP5G_7dpf_LongStim_F1_RESULTS_Selma.mat']);

n_trial=1; %trial. 3 in total
n_cells=size(Fish{1,1}.DSdeltaFoF{n_trial},1);
K=4 %stimulus repetition
theta=[pi/2 ; 0 ; -pi/2 ; -pi];

DS_DeltaFoF=cell(1,4); %une cellule par angle theta
DS_DeltaFoFtotal=zeros(n_cells,4); %for each i=1:4 mean of the response to theta_i oriented moving bar stimulus
for j=1:4

for n=1:n_cells
    for k=1:K
DS_DeltaFoF{j}(n,k)=Fish{1,1}.DSdeltaFoF{n_trial}{n,k}(j);
DS_DeltaFoFtotal(n,j)=Fish{1,1}.DSdeltaFoFtotal{n_trial}{1,n}(j);
    end
    myvar(n,j)=var(DS_DeltaFoF{j}(n,:));
end
end
%% === statistic on direction response ==================================%%
close all
for n=1:n_cells
R(n)=sum(abs(DS_DeltaFoFtotal(n,:)-min(min(DS_DeltaFoFtotal(:,:)))).*exp(1i*theta(:))')/sum(abs(DS_DeltaFoFtotal(n,:)-min(min(DS_DeltaFoFtotal(:,:)))'));
normR(n)=sqrt(real(R(n))^2+imag(R(n))^2);
V(n)=1-norm(R(n));

circular_tc(n,:)=abs(DS_DeltaFoFtotal(n,:)-min(min(DS_DeltaFoFtotal(:,:)))).*exp(1i*theta(:))'/sum(abs(DS_DeltaFoFtotal(n,:)-min(min(DS_DeltaFoFtotal(:,:)))'));
mymean(n)=mean(circular_tc(n,:)); %autre moyen de calculer le R


Rsym(n)=sum((DS_DeltaFoFtotal(n,:)-min(min(DS_DeltaFoFtotal(:,:)))).*exp(1i*mod(theta(:),pi))')/sum(abs(DS_DeltaFoFtotal(n,:)-min(min(DS_DeltaFoFtotal(:,:)))'));
normRsym(n)=sqrt(real(Rsym(n))^2+imag(Rsym(n))^2);
Vsym(n)=1-norm(Rsym(n));
end

%% ==== plot circular tuning curves
figure;
hold on;
mycolors=['g','b','v','c','y','r'];
for n=1:n_cells
plot((1-V(n))*real(circular_tc(n,:)),(1-V(n))*imag(circular_tc(n,:)),'x');%,mycolors(mod(n,6)+1)); 
end
%% ==== plot angular responses
figure;
plot(real(R(:)),imag(R(:)),'x'); title('orientation selectivity');
figure; 
hold on
for n=1:n_cells 
    plot(real(mymean(n)),imag(mymean(n)),'x') 
end
hold off
title('orientation selectivity');
figure;
plot(real(Rsym(:)),imag(Rsym(:)),'x'); title('direction selectivity');
figure; 
hold on
for n=1:n_cells 
    plot(real(Rsym(n)),imag(Rsym(n)),'x') 
end
hold off
title('direction selectivity');
figure;
plot(V);
hold on; 
plot(mean(V)*ones(size(V)),'c');
hold off;
title('circular variance for orientation selectivity')
figure;
plot(Vsym)
title('circular variance for direction selectivity')

%% === plot taking into account neuron selectivity
close all

theta_temoin=linspace(-pi,pi,105);

figure; 
plot(cos(theta_temoin),sin(theta_temoin),'k')
hold on
for n=1:n_cells 
    plot(V(n)*1/norm(mymean(n))*real(mymean(n)),V(n)*1/norm(mymean(n))*imag(mymean(n)),'x') 
end
hold off
title('orientation selectivity modulated by circular variance');

figure; 
plot(cos(theta_temoin),sin(theta_temoin),'k')
hold on
for n=1:n_cells 
    plot(Vsym(n)*1/norm(Rsym(n))*real(Rsym(n)),Vsym(n)*1/norm(Rsym(n))*imag(Rsym(n)),'x') 
end
hold off
title('direction selectivity modulated by circular variance');

%% find prefered directions


DS=sum(real(R'*exp(-1i*theta_temoin)));
figure;
plot(real(DS(:)'.*exp(1i*theta_temoin)),imag(DS(:)'.*exp(1i*theta_temoin)),'x')
% hold on; 
% plot(real(R(:)),imag(R(:)),'x'); 
% plot(real(DS(1:26:105).*exp(1i*theta_temoin(1:26:105))),imag(DS(1:26:105).*exp(1i*theta_temoin(1:26:105))),'o')
% hold off;


%% === TEMPLATE MATCHING
n_angles=size(theta,1);
n_stim_total=K*n_angles;
template=zeros(n_angles,n_cells);
actual=DS_DeltaFoF;
% DS_DeltaFoF=cell(1,4); %%%%%%%%%--->>> une cellule par angle theta
Similarity=zeros(n_angles,n_stim_total);


for j=1:K
    for n=1:n_cells
template(j,n)=DS_DeltaFoFtotal(n,j)';
    end
end

for nangle=1:n_angles
for k=1:K
Similarity(nangle,k*nangle)=sum(template(nangle,:)'.*actual{k}(:,nangle))/(norm(template(nangle,:))*norm(actual{k}(:,nangle)));
end
end

for nangle=1:n_angles
figure;
hold on;
for j=1:n_stim_total
plot(real(Similarity(nangle,j)*exp(1i*theta(nangle))),imag(Similarity(nangle,j)*exp(1i*theta(nangle))),'o')
end
title(['Similarity for angle ',num2str(nangle),'.']);
hold off
end

%% === map des sensibilités =============================================%%

myangles=angle(R);%mymean);%atan(imag(mymean)./real(mymean));

% figure;
hist(myangles);%atan(imag(mymean)./real(mymean)));
title('histogramme des angles modulo pi');

pos=Fish{1,1}.AllCells{1,1}.cell;
d2pos=cell(n_cells,1);
for n=1:n_cells
y=floor(pos{n}/249)+1;
x=mod(pos{n},249)+1;
for i=1:numel(y)
if x(i)==0
x=249;
end
end
d2pos{n}=[x,y];
end

mybkg=Fish{1,1}.Raster{1,1}.imageAvg;
bkg=zeros(size(mybkg,1),size(mybkg,2),3);
for x=1:size(mybkg,2)
    for y=1:size(mybkg,1)
bkg(y,x,:)=[0,0,0.01*mybkg(y,x)];
    end
end

Bin=linspace(min(myangles),max(myangles),100);
color=jet(numel(Bin));

figure;

subplot(5,1,1:4)
imshow(rescalegd(mybkg));

mymap=0.5*ones(size(bkg,1),size(bkg,2),3);%zeros(size(bkg));%
 A=zeros(size(bkg,1),size(bkg,2));
for n=1:n_cells
for j=1:size(pos{n},1)
[mymin,myidx]=min(abs(Bin-myangles(n)*ones(size(Bin))));
mymap(d2pos{n}(j,1),d2pos{n}(j,2),:)=color(myidx,:);
end
end
A(Fish{1,1}.AllCells{1,1}.bkg>0.0) = 1;
hold on;
hshow=imshow(mymap);
set(hshow,'AlphaData',A)
hold off
subplot(5,1,5)
mycolorbar=zeros(10,size(color,1),size(color,2));
for j=1:10
mycolorbar(j,:,:)=color;
end
imshow(mycolorbar)
title('Prefered direction U-R-D-L');


%% === map des variances=================================================%%

for j=1:4
    


figure;
hist(myvar(:,j));%atan(imag(mymean)./real(mymean)));
title(['distribution de variance pour stimulus ', num2str(j),'.']);

pos=Fish{1,1}.AllCells{1,1}.cell;
d2pos=cell(n_cells,1);
for n=1:n_cells
y=floor(pos{n}/249)+1;
x=mod(pos{n},249)+1;
for i=1:numel(y)
if x(i)==0
x=249;
end
end
d2pos{n}=[x,y];
end




mybkg=Fish{1,1}.Raster{1,1}.imageAvg;
bkg=zeros(size(mybkg,1),size(mybkg,2),3);
for x=1:size(mybkg,2)
    for y=1:size(mybkg,1)
bkg(y,x,:)=[0,0,0.01*mybkg(y,x)];
    end
end

Bin=linspace(min(myvar(:,j)),max(myvar(:,j)),100);
color=jet(numel(Bin));

figure;

subplot(5,1,1:4)
imshow(rescalegd(mybkg));

mymap=0.5*ones(size(bkg,1),size(bkg,2),3);%zeros(size(bkg));%
 A=zeros(size(bkg,1),size(bkg,2));
for n=1:n_cells
for l=1:size(pos{n},1)
[mymin,myidx]=min(abs(Bin-myvar(n,j)*ones(size(Bin))));
mymap(d2pos{n}(l,1),d2pos{n}(l,2),:)=color(myidx,:);
end
end
A(Fish{1,1}.AllCells{1,1}.bkg>0.0) = 1;
hold on;
hshow=imshow(mymap);
set(hshow,'AlphaData',A)
hold off
subplot(5,1,5)
mycolorbar=zeros(10,size(color,1),size(color,2));
for l=1:10
mycolorbar(l,:,:)=color;
end
imshow(mycolorbar)
title(['Variance per neuron for stimulus ',num2str(j),'.']);


end

%% === big map =========================================================%%

for j=1:4
    

pos=Fish{1,1}.AllCells{1,1}.cell;
d2pos=cell(n_cells,1);
for n=1:n_cells
y=floor(pos{n}/249)+1;
x=mod(pos{n},249)+1;
for i=1:numel(y)
if x(i)==0
x=249;
end
end
d2pos{n}=[x,y];
end




mybkg=Fish{1,1}.Raster{1,1}.imageAvg;
bkg=zeros(size(mybkg,1),size(mybkg,2),3);
for x=1:size(mybkg,2)
    for y=1:size(mybkg,1)
bkg(y,x,:)=[0,0,0.01*mybkg(y,x)];
    end
end

Bin=linspace(min(myvar(:,j)),max(myvar(:,j)),100);
color=jet(numel(Bin));

figure;

subplot(5,1,1:4)
imshow(rescalegd(mybkg));

mymap=0.5*ones(size(bkg,1),size(bkg,2),3);%zeros(size(bkg));%
 A=zeros(size(bkg,1),size(bkg,2));
for n=1:n_cells
for l=1:size(pos{n},1)
[mymin,myidx]=min(abs(Bin-myvar(n,j)*ones(size(Bin))));
mymap(d2pos{n}(l,1),d2pos{n}(l,2),:)=color(myidx,:);
end
end
A(Fish{1,1}.AllCells{1,1}.bkg>0.0) = 1;
hold on;
hshow=imshow(mymap);
set(hshow,'AlphaData',A)
hold off
subplot(5,1,5)
mycolorbar=zeros(10,size(color,1),size(color,2));
for l=1:10
mycolorbar(l,:,:)=color;
end
imshow(mycolorbar)
title(['Variance per neuron for stimulus ',num2str(j),'.']);


end
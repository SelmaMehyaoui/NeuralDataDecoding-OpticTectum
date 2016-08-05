%%===========================decoding optic tectum rf====================%%

clear all; 
close all; 
clc
mypath='/users/zfne/mehyaoui/Desktop/phD_Lab2/Donnees Optic Tectum/';
load([mypath,'141024_GCaMP5G_7dpf_LongStim_F1/Results/141024_GCaMP5G_7dpf_LongStim_F1_RESULTS_Selma.mat']);

n_trial=1; %trial. 3 in total
n_cells=size(Fish{1,1}.RFdeltaFoF{n_trial},1);
K=10; %stimulus repetition (en fait 20 mais deux fois dix tirages dc on moyenne sur stimulus identiques)
n_stim=720;  %nb de stimulus (chaque stim 

RF_DeltaFoF=Fish{1,1}.RFdeltaFoF{n_trial}; %une cellule par stimulation et par neurone, chaque cellule taille 9x4
RF_DeltaFoFtotal=Fish{1,1}.RFdeltaFoFtotal{n_trial};
Stimulus=Fish{1,1}.Timing{1,1}.randomOrders;

TotalMeanResponse=zeros(n_cells,36);
M=cell(n_cells,36); %M : ReOrganizedResponses each cell contains 20 dim vector 

for n=1:n_cells
for k=1:36
for j=1:20

new_j=mod(j,10);
if j==10
new_j=10;
elseif j==20
new_j=10;
end

% [~,idx]=find(Fish{1,1}.Timing{1,1}.randomOrders{1,new_j}==k);
idx=Fish{1,1}.Timing{1,1}.randomOrders{1,new_j}(k);
M{n,k}(j)=RF_DeltaFoF{n,j}(idx);
end
myvar(n,k)=var(M{n,k});
end
end

for n=1:n_cells
my_rf=RF_DeltaFoFtotal{n}';
TotalMeanResponse(n,:)=my_rf(:)';
end

rpsformax=TotalMeanResponse';

[mymax,idmax]=max(rpsformax);

%% === plot rough receptive fields ===================================%%

close all

imagesc((TotalMeanResponse-min(min(TotalMeanResponse))).^2);
title('Responses neuron per neuron');

mystimulus=sort(Stimulus{1});
figure;
hold on;
for n=1:n_cells
plot(mystimulus,TotalMeanResponse(n,:));

end
title('Tuning curves');
hold off

figure; imagesc(cov((TotalMeanResponse-min(min(TotalMeanResponse))*ones(size(TotalMeanResponse))))');
% figure; 
% hold on; 
% for n=1:n_cells
% f = fit(mystimulus.',TotalMeanResponse(n,:).','gauss2');
% plot(f,mystimulus,TotalMeanResponse(n,:))
% end
% hold off;

%% === map des RF ======================================================%%

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

Bin=linspace(min(idmax),max(idmax),100);
color=jet(numel(Bin));


figure;
imshow(rescalegd(mybkg));

mymap=0.5*ones(size(bkg,1),size(bkg,2),3);%zeros(size(bkg));%
 A=zeros(size(bkg,1),size(bkg,2));
for n=1:n_cells
for j=1:size(pos{n},1)
mymap(d2pos{n}(j,1),d2pos{n}(j,2),:)=color(idmax(n),:);
end
end
A(Fish{1,1}.AllCells{1,1}.bkg>0.0) = 1;
hold on;
hshow=imshow(mymap);
set(hshow,'AlphaData',A)
hold off


%% === map des variances de RF===========================================%%
for j=1:36
    
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
%% === statistic on RF =================================================%%
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


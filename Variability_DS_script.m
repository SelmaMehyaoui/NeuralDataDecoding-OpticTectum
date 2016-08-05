%% ===================Optic Tectum : Variability in Direction Decoding == %%
% Moving bars. Four cardinal directions

clear all;
close all; 
clc
mypath='/users/zfne/mehyaoui/Desktop/phD_Lab2/Donnees Optic Tectum/';
load([mypath,'141024_GCaMP5G_7dpf_LongStim_F1/Results/141024_GCaMP5G_7dpf_LongStim_F1_RESULTS_Selma.mat']);
n_cells=312;

%% ================== Defining variables =============================== %%
n_trial=1;
%n_cells=size(Fish{1,1}.RFdeltaFoF{n_trial},1);
DS=zeros(n_cells, 1);
n_stim=4; %nb de répétitions d'un stimulus donné
N_orientations=4; %URDL
N_stim=N_orientations*n_stim; %nb total de stimulations


DS_DeltaFoF=cell(1,N_orientations); %une cellule par angle theta
DS_DeltaFoFtotal=zeros(n_cells,N_orientations); %for each i=1:4 mean of the response to theta_i oriented moving bar stimulus
for j=1:N_orientations
for n=1:n_cells
    for k=1:n_stim
DS_DeltaFoF{j}(n,k)=Fish{1,1}.DSdeltaFoF{n_trial}{n,k}(j);
DS_DeltaFoFtotal(n,j)=Fish{1,1}.DSdeltaFoFtotal{n_trial}{1,n}(j);
    end
    myvar(n,j)=var(DS_DeltaFoF{j}(n,:));
end
end

%% ================== Noise correlations =============================== %%
Noise_Corr=cell(n_stim,1);
for k=1:n_stim
Noise_Corr{k}=corrcoef(DS_DeltaFoF{k}');
my_Cov=cov(DS_DeltaFoF{k}');
Noise_Corr{k}(isnan(Noise_Corr{k}))=2*sign(my_Cov(isnan(Noise_Corr{k})));
end

figure;
for k=1:n_stim
subplot(n_stim,1,k)
imagesc(Noise_Corr{k})
colorbar
title(['Noise correlation for stimulus #',num2str(k),'.']);
end

meanNC=zeros(size(Noise_Corr{1}));
for k=1:n_stim
meanNC=meanNC+Noise_Corr{k};
end
meanNC=1/4*meanNC;
% my_Cov=cov(DS_DeltaFoF{k}');
% meanNC(isnan(meanNC))=2*sign(my_Cov(isnan(meanNC)));

figure;
imagesc(meanNC)
colorbar
title('Noise Correlations')

Znoise = linkage(meanNC);
c=1.153
Tnoise = cluster(Znoise,'cutoff',c);
colorbar;

nNC_clust=max(Tnoise)

[Tsort,Idx]=sort(Tnoise);

Sort_meanNC=meanNC(Idx,Idx);
figure; imagesc(Sort_meanNC); colorbar;
title('Sorted Noise correlations');

%% ================== Signal correlations ============================== %%

Signal_Corr=corrcoef(DS_DeltaFoFtotal');

figure;
imagesc(Signal_Corr)
colorbar
title('Signal Correlations')

Zsignal = linkage(Signal_Corr);
c=0.12;
Tsignal = cluster(Zsignal,'cutoff',c);
colorbar;

nSC_clust=max(Tsignal)

[Tsort,Idx]=sort(Tsignal);

Sort_Signal_Corr=Signal_Corr(Idx,Idx);
figure; imagesc(Sort_Signal_Corr); colorbar;
title('Sorted Signal correlations');
%% ================= 
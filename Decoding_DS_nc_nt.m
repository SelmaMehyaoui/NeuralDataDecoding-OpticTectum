%% ===================Optic Tectum : Orientation decoding ============== %%
function [Linear_Decoder_Error, Linear_Decoder_Precision, Match_Decoder_Error,Precision4_Decoder_Error,Precision_Decoder_Error,Match4_Decoder_Error,Linear_Decoder2_Error, Mean4_Decoder2_Error, Pop_Vector_Error,Template_Matching_Error,ZS_Template_Matching_Error,ML_Error]=Decoding_DS_nc_nt(id_cells, id_train, id_test, ifplots,whatdecoder,saves_struct_number,tol1,tol2)
% Moving bars. Four cardinal directions

% id cells, id train, id tests are binary vectors

% whatdecoder : binary vector indicated what decoder are built and tested
% whatdecoder(1) decodeur lineaire (via Moore Penrose pseudo-inverse)
% whatdecoder(2) decodeur lineaire sur le vecteur (cos x sinx cos2x sin2x)
% whatdecoder(3) decodeur lineaire2 (en argmax simple sur les distributions)
% whatdecoder(4) decodeur lineaire2 sur le vecteur (cosx sinx cos2x sin2x)
% whatdecoder(5) population vector decoding
% whatdecoder(6) template matching
% whatdecoder(7) normalized template matching
% whatdecoder(8) bayesian maximum likelihood 

% si whatdecoder(i) est ? zero le decodeur n'est pas construit et un
% vecteur d'erreur r?duit au singleton z?ro est renvoy?

% si whatdecoder(i) est ? un le decodeur est construit et le vecteur
% d'erreur est renvoy? en sortie de la fonction

% clear all;
% 
% close all; 
mypath='/Users/selmamehyaoui/Desktop/Optic Tectum Data/Donnees John/';
load([mypath,'141024_GCaMP5G_7dpf_LongStim_F1/Results/141024_GCaMP5G_7dpf_LongStim_F1_RESULTS_Selma.mat']);

%% ================== Set default output to zero ======================= %%

Linear_Decoder_Error=0;
Match4_Decoder_Error=0;
Linear_Decoder2_Error=0;
Mean4_Decoder2_Error=0;
Pop_Vector_Error=0;
Template_Matching_Error=0;
ZS_Template_Matching_Error=0;
ML_Error=0;

%% ================== Defining variables =============================== %%
n_trial=1;
n_cells_tot=size(Fish{1,1}.RFdeltaFoF{n_trial},1)
n_cells=sum(id_cells);

n_stim=4; %nb de repetitions d'un stimulus donn??
N_orientations=4; %URDL % sin2theta ueless ici !! 
N_stim=N_orientations*n_stim; %nb total de stimulations

n_train=sum(id_train); 
n_test=sum(id_test); 

DS=zeros(n_cells, 1);

DS_DeltaFoF=cell(1,4); 
%une cellule par angle theta

DS_DeltaFoFtotal=zeros(n_cells,4); 
%for each i=1:4 mean of the response to theta_i oriented moving bar stimulus

Ordered_DS_DeltaFoF=zeros(n_cells,16);
DS_DeltaFoFtotaltrain=zeros(n_cells,4);

for j=1:N_orientations

size(id_cells)
size((1:n_cells_tot))
myneworder=id_cells.*(1:n_cells_tot);
myneworder(myneworder==0)=[];

for n=1:n_cells
    for k=1:n_stim
DS_DeltaFoF{j}(n,k)=Fish{1,1}.DSdeltaFoF{n_trial}{myneworder(n),k}(j);
DS_DeltaFoFtotal(n,j)=Fish{1,1}.DSdeltaFoFtotal{n_trial}{1,myneworder(n)}(j);
    end
    myvar(n,j)=var(DS_DeltaFoF{j}(n,:));
end
end


for j=1:N_orientations
for k=1:n_stim
for n=1:n_cells
Ordered_DS_DeltaFoF(n,(j-1)*4+k)=DS_DeltaFoF{j}(n,k);
end
end
end

for n=1:n_cells
a(n,:)=Ordered_DS_DeltaFoF(n,:).*id_train;
end
a(:,id_train==0)=[];
DS_DeltaFoF_train=a;

for j=1:N_orientations
temoinfor{j}=zeros(1,16);
temoinfor{j}((j-1)*4+1:j*4)=1;
temoinfor{j}(id_train==0)=[];
end



for j=1:N_orientations
DS_DeltaFoFtotaltrain(n,j)=mean(a(n,:).*temoinfor{j});
end

%defining DeltaFof for tested trials
for n=1:n_cells
A(n,:)=Ordered_DS_DeltaFoF(n,:).*id_test;
end
A(:,id_test==0)=[];


%idtrain idtest
% 4 repetitions de chacun des 4 stimuli soit 16 trials au total
my_stim=[ones(1,4),2*ones(1,4),3*ones(1,4),4*ones(1,4)];
my_stim_train=my_stim.*id_train;
my_stim_train(my_stim_train==0)=[];
my_stim=my_stim.*id_test
my_stim(my_stim==0)=[];

cos_stim_train=round(cos(pi/2-(my_stim_train-1)*pi/2));
sin_stim_train=round(sin(pi/2-(my_stim_train-1)*pi/2));
cos_2stim_train=round(cos(pi/2-(my_stim_train-1)*pi/2));
sin_2stim_train=round(sin(pi/2-(my_stim_train-1)*pi/2));
cossincos2sin2_train=[cos_stim_train;sin_stim_train;cos_2stim_train;sin_2stim_train];

cos_stim=round(cos(pi/2-(my_stim-1)*pi/2));
sin_stim=round(sin(pi/2-(my_stim-1)*pi/2));
cos_2stim=round(cos(pi/2-(my_stim-1)*pi/2));
sin_2stim=round(sin(pi/2-(my_stim-1)*pi/2));
cossincos2sin2_stim=[cos_stim;sin_stim;cos_2stim;sin_2stim]


%% ======================= Decodeur lineaire ============================%%

% whatdecoder(1) decodeur lineaire (via Moore Penrose pseudo-inverse)

if whatdecoder(1)
W=my_stim_train*pinv(DS_DeltaFoF_train,tol1);
Decoded_linear=W*A;

Linear_Decoder_Error=sum(1-round(Decoded_linear)==my_stim)/n_test*100;
Linear_Decoder_Precision=sum((Decoded_linear-my_stim).^2)/n_test;
end

%% ============= Decodeur lin?aire en cos sin cos2 sin2==================%%

% whatdecoder(2) decodeur lineaire
if whatdecoder(2)
myinv=pinv(DS_DeltaFoF_train,tol2);
W4=cossincos2sin2_train*myinv;
Decoded_linear4=W4*A;
Decoded_linear4=Decoded_linear4./max(max(Decoded_linear4));
for j=1:n_test
f = @(x) norm([cos(x);sin(x);cos(2*x);sin(2*x)]-Decoded_linear4(:,j))^2;
options = optimoptions('fminunc','Algorithm','quasi-newton');
x0=acos(Decoded_linear4(1,j));
[x(j), fval, exitflag, output] = fminunc(f,x0,options);
end

Precision_Decoder_Error=(2*(1-x/pi)-my_stim).^2
Match_Decoder_Error=1-(round(2*(1-x/pi))==my_stim);

Precision4_Decoder_Error=sum((Decoded_linear4-cossincos2sin2_stim).^2);
Match4_Decoder_Error=(sum(Decoded_linear4==cossincos2sin2_stim)<size(Decoded_linear4,1));

Match_Decoder_Error=sum(Match_Decoder_Error)/n_test*100;
Precision4_Decoder_Error=sum(Precision4_Decoder_Error)/n_test;
Precision_Decoder_Error=sum(Precision_Decoder_Error)/n_test;
Match4_Decoder_Error=sum(Match4_Decoder_Error)/n_test*100; 

end

%% ======================= Decodeur lin?aire2 ===========================%%

% whatdecoder(3) decodeur lineaire

if whatdecoder(3)

[~,PrefDir]=max((DS_DeltaFoFtotaltrain)');%./myvar)');
Linear_Decoder2_Error=0;

for l=1:n_test

Linearly_decoded2_Dir(l)=sum(A(:,l).*(pi/2-(PrefDir-1)*pi/2)')/sum(A(:,l));
Linear_Decoder2_Error=Linear_Decoder2_Error+1-isequal(round(-(Linearly_decoded2_Dir(l)-pi/2)/pi*2+1),my_stim(l));
   
end

Linear_Decoder2_Error=Linear_Decoder2_Error/n_test*100; %error in percentage

if ifplots
figure; hold on; plot(my_stim,'c*','linewidth',2); plot(round(-(Linearly_decoded2_Dir-pi/2)/pi*2+1),'bx');  hold off
disp(['Error = ',num2str(Pop_Vector_Error),' %']);
title('Decoding with Linear decoder 2');
xlabel('trial')
ylabel('decoded stimulus (blue) , actual stimulus (red) ')
end

end

%% ============ Decoding lineaire2 en cos sin cos2 sin2 ================ %%

% whatdecoder(4) decodeur lineaire
if whatdecoder(4)
Prefcos=zeros(3,n_cells) ; %dans l'experience consid?r?e cos =-1 0 1
Prefsin=zeros(3,n_cells); %dans l'experience consid?r?e sin=-1 0 1
Prefcos2=zeros(2,n_cells); %dans l'exp?rience consid?r?e cos2= -1 ou 1
%Prefsin2= xxx ; %peu ns importe ici, sin2=0 pr tous les stimuli !

whatjmeans=[-1 0 1];
whatkmeans=[-1 1];

for j=1:3
for n=1:n_cells
Prefcos(j,n)=sum(Ordered_DS_DeltaFoF(n,cos_stim==whatjmeans(j)))/sum(cos_stim==whatjmeans(j));
Prefsin(j,n)=sum(Ordered_DS_DeltaFoF(n,sin_stim==whatjmeans(j)))/sum(sin_stim==whatjmeans(j));
end
end
[~,MyPrefcos]=max((Prefcos));
MyPrefcos(MyPrefcos==1)=-1;
MyPrefcos(MyPrefcos==2)=0;
MyPrefcos(MyPrefcos==3)=1;

[~,MyPrefsin]=max((Prefsin));
MyPrefsin(MyPrefsin==1)=-1;
MyPrefsin(MyPrefsin==2)=0;
MyPrefsin(MyPrefsin==3)=1;
for k=1:2
for n=1:n_cells
Prefcos2(k,n)=sum(Ordered_DS_DeltaFoF(n,cos_2stim==whatkmeans(k)))/sum(cos_2stim==whatkmeans(k));
end
end
[~,MyPrefcos2]=max((Prefcos2));
MyPrefcos2(MyPrefcos2==1)=-1;
MyPrefcos2(MyPrefcos2==2)=1;
MyPrefsin2=zeros(size(MyPrefcos2));
MyPrefcossin12=[MyPrefcos;MyPrefsin;MyPrefcos2;MyPrefsin2];

Linear4_Decoder2_Error=zeros(4,1);

for l=1:n_test
for k=1:4
Linearly4_decoded2_Dir(l,k)=sum(A(:,l).*MyPrefcossin12(k,:)')/sum(A(:,l));
Linear4_Decoder2_Error(k)=Linear4_Decoder2_Error(k)+1-isequal(round(Linearly4_decoded2_Dir(l,k)),cos_stim(l));
end
end

Linear4_Decoder2_Error=Linear4_Decoder2_Error./n_test*100; %error in percentage
Mean4_Decoder2_Error=mean(Linear4_Decoder2_Error);

end

%% ================== Population Vector Decoding ======================= %%
% whatdecoder(5) population vector
if whatdecoder(5)

[~,PrefDir]=max((DS_DeltaFoFtotaltrain)');%./myvar)');
 
Pop_Vector_Error=0;

%Linear4_Decoder_Error ds la boucle for qui vient

for l=1:n_test
decoded_Dir2(l)=angle(sum(A(:,l).*exp(i*(pi/2-(PrefDir-1)*pi/2))')/sum(A(:,l)));
Pop_Vector_Error=Pop_Vector_Error+1-isequal(round(-(decoded_Dir2(l)-pi/2)/pi*2+1),my_stim(l));
end
Pop_Vector_Error=Pop_Vector_Error/n_test*100; %error in percentage

if ifplots

figure; hold on; plot(my_stim,'c*','linewidth',2); plot(round(-(decoded_Dir2-pi/2)/pi*2+1),'bx');  hold off
disp(['Error = ',num2str(Pop_Vector_Error),' %']);
title('Decoding with Population Vector');
xlabel('trial')
ylabel('decoded stimulus (blue) , actual stimulus (red) ')
end

end

%% ================== Template Matching Decoding ======================== %%
% whatdecoder(6) template matching

if whatdecoder(6)

Response_Template=DS_DeltaFoFtotaltrain;


Template_Matching_Error=0;

%here change the number of trials
for l=1:n_test
for theta=1:4
decoded_Dir_TM(l,theta)=sum(A(:,l).*Response_Template(:,theta))/(norm(A(:,l)))*norm(Response_Template(:,theta));
end
end

[~,Decoded_Dir_TM]=max(decoded_Dir_TM');

% here change the number of trials
for l=1:n_test
Template_Matching_Error=Template_Matching_Error+1-isequal(Decoded_Dir_TM(l),my_stim(l)); %error in percentage
end
Template_Matching_Error=Template_Matching_Error*100/n_test;

if ifplots

disp(['Error =',num2str(Template_Matching_Error),' %'])
figure; hold on; plot(my_stim,'c*','linewidth',2); plot(Decoded_Dir_TM,'bx'); hold off
title('Decoding with Template Matching');
xlabel('trial')
ylabel('decoded stimulus (blue) , actual stimulus (red) ')

end

end

%% =========== Normalized Pattern Matching Decoding ==================== %%
% whatdecoder(7) normalized template matching

%The only difference is that vectors are z-scored

if whatdecoder(7)

ZS_Response_Template=zscore(DS_DeltaFoFtotaltrain')';
ZS_Ordered_DS_DeltaFoF=zscore(A')';

ZS_Template_Matching_Error=0;

%here change number of trials
for l=1:n_test
for theta=1:4
ZS_decoded_Dir_TM(l,theta)=sum(ZS_Ordered_DS_DeltaFoF(:,l).*ZS_Response_Template(:,theta))/(norm(ZS_Ordered_DS_DeltaFoF(:,l)))*norm(ZS_Response_Template(:,theta));
end
end

[~,ZS_Decoded_Dir_TM]=max(ZS_decoded_Dir_TM');

%here change the number of trials
for l=1:n_test
ZS_Template_Matching_Error=ZS_Template_Matching_Error+1-isequal(ZS_Decoded_Dir_TM(l),my_stim(l)); %error in percentage
end
ZS_Template_Matching_Error=ZS_Template_Matching_Error*100/n_test;

if ifplots
disp(['Error =',num2str(ZS_Template_Matching_Error),' %'])

figure; hold on; plot(my_stim,'c*','linewidth',2); plot(ZS_Decoded_Dir_TM,'bx');  hold off
title('Decoding with Normalized Template Matching')
xlabel('trial')
ylabel('decoded stimulus (blue) , actual stimulus (red) ')

end

end

%% ============= Bayesian Maximum-Likelihood =========================== %%

% whatdecoder(8) bayesian maximum likelihood 

% P(theta|Apop)=[P(Apop|theta).P(theta)]/P(Apop)
% %posterior=likelihood.prior/prob(act pattern)
% %P(theta|Apop) ~ Produit_[i=1:n](Ptheta|A_i)
% %P(A_i|theta)P(theta)/P(A_i)
% * fit gaussien sur distri de Ai ?? theta fix??
% ** p(theta) ?? 1/4
% *** P(Ai) = n_occ(Ai)/16


% "ONE STIMULUS PROB is estimated from training set trials where THIS stimulus direction was
% presented"

% idea : enregistrer dans une variable valise le nombre de fois o? le
% stimulus a ?t? pr?sent? / nbtotalpresprunstim. Compris entre0 et 1.Plus
% on est proche de zero plus il est surprenantd'avoir des r?sultats, plus 
% on est proche de 1, plus un bon score est attendu.
if whatdecoder(8)

for j=1:4 % j : theta
P_theta(j)=1/4;    
for n=1:n_cells
[y,x]=hist(Ordered_DS_DeltaFoF(n,(j-1)*4+1:j*4).*id_train((j-1)*4+1:j*4),100);
y=y/sum(y);
for k=1:n_test % repet indpendemment du stimulus en question
[~,myx]=min(abs(A(n,k)-x));
epsilon=0.0000001; %10-7
if y(myx)==0
P_A_theta{j}(k,n)=epsilon;
else
P_A_theta{j}(k,n)=y(myx);
end
P(n,k)=sum(DS_DeltaFoF_train(n,:)==A(n,k))/n_train;

P_theta_pop_n{j}(n,k)=P_A_theta{j}(k,n)*P_theta(j)/P(n,k);
end
end

%here change the number of trials
for k=1:n_test
P_theta_pop(j,k)=prod(P_theta_pop_n{j}(:,k));
end
end


[~,ML_decoded_Dir]=max(log(P_theta_pop));

ML_Error=(n_test-sum((my_stim-ML_decoded_Dir)==0))*100/n_test; % en pourcentage

if ifplots
disp(['Error = ',num2str(ML_Error),' %'])
figure; hold on; plot(my_stim,'c*','linewidth',2); plot(ML_decoded_Dir,'bx');  hold off
title('Decoding with Bayesian ML');
xlabel('trial')
ylabel('decoded stimulus (blue) , actual stimulus (red) ')
end

if saves_struct_number>0
save(['Mydecoders_',num2str(saves_struct_number)],'MyPrefcossin12','PrefDir','Response_Template','ZS_Response_Template','P_theta_pop_n');
end
end
%% =====================THIS IS THE END ================================ %%

end
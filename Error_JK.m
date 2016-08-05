%% ATTENTION, SCRIPT NEURONE n SPECIFIQUE, ie NEURONE n FIXE
%% D AILLEURS A TRANSFORMER EN FONCTION DE n ULTERIEUREMENT

%% SCRIPT FOR COMPUTING NEURON I DECODING IMPROVEMENT
%ici n_cells renvoie a la taille de la sous pop a laquelle l apport du neurone i est calcule
 
function [Npopconfig,Distance,mylineardecoder_error,mylinear4decoder_error,mypopvector_error,myTM_error,myNTM_error,myBML_error,JK_index_linear,JK_index_linear4,JK_index_popvector,JK_index_TM,JK_index_NTM,JK_index_BML]=Error_JK(My_n, ifplots)
%My_n = neurone dont on p?se l'apport aux sous-pop dont il fait
%potentiellement partie

Npopconfig=zeros(1,312-1);

for n_cells=2:1:312

Npopconfig(n_cells)=nchoosek(312,n_cells-1);
myNpopconfig=Npopconfig(n_cells);

for popconfig=1:myNpopconfig

c_cells = cvpartition(312,'HoldOut',n_cells-1);
%on fait le pari qu'? chaque tirage ?chantillon diff?rent

%{popconfig}

id_pop=[zeros(1,My_n-1),1,zeros(1,312-My_n)];
id_cells=c_cells.test'+id_pop;
id_cells_bis=id_cells;
id_cells_bis(My_n)=0;
%? changer, cette fois la composition de la population importe
%elle doit contenir le neurone i

%calcul de la distance du neurone i au reste de la sous pop

d1=0;
compteur=0;
while (compteur==0)&&(My_n-d1>0)
d1=d1+1;
if d1==My_n
d1=400;
break
end

idforc=c_cells.test';
compteur=sum(idforc(My_n-d1:My_n));
end
d1=d1-1;

d2=0;
compteur=0;
while (compteur==0)&&(My_n+d2<312)
d2=d2+1;
idforc=c_cells.test';
compteur=sum(idforc(My_n:My_n+d2));
end
d2=d2-1;


Distance(n_cells,popconfig)=min(d1,d2);

for p=0.1:0.1:1
k=int8(p*10);
c=cvpartition(16,'Holdout',1/16);%p mais 1/16 pour un one left out
id_test=test(c);
id_train=training(c);
ifplots=0;
[MyLinear_Decoder_error(n_cells,k),MyLinear4_Decoder_error(n_cells,k),Mypopvector_error(n_cells,k),MyTM_error(n_cells,k),MyNTM_error(n_cells,k),MyBML_error(n_cells,k)]=Decoding_DS_nc_nt(id_cells, id_train', id_test', ifplots);
[MyLinear_Decoder_error_bis(n_cells,k),MyLinear4_Decoder_error_bis(n_cells,k),Mypopvector_error_bis(n_cells,k),MyTM_error_bis(n_cells,k),MyNTM_error_bis(n_cells,k),MyBML_error_bis(n_cells,k)]=Decoding_DS_nc_nt(id_cells_bis, id_train', id_test', ifplots);
end
mylineardecoder_error(n_cells,popconfig)=mean(MyLinear_Decoder_error(n_cells,:));
mylinear4decoder_error(n_cells,popconfig)=mean(MyLinear4_Decoder_error(n_cells,:));
mypopvector_error(n_cells,popconfig)=mean(Mypopvector_error(n_cells,:));
myTM_error(n_cells,popconfig)=mean(MyTM_error(n_cells,:));
myNTM_error(n_cells,popconfig)=mean(MyNTM_error(n_cells,:));
myBML_error(n_cells,popconfig)=mean(MyBML_error(n_cells,:));

mylineardecoder_error_bis(n_cells,popconfig)=mean(MyLinear_Decoder_error_bis(n_cells,:));
mylinear4decoder_error_bis(n_cells,popconfig)=mean(MyLinear4_Decoder_error_bis(n_cells,:));
mypopvector_error_bis(n_cells,popconfig)=mean(Mypopvector_error_bis(n_cells,:));
myTM_error_bis(n_cells,popconfig)=mean(MyTM_error_bis(n_cells,:));
myNTM_error_bis(n_cells,popconfig)=mean(MyNTM_error_bis(n_cells,:));
myBML_error_bis(n_cells,popconfig)=mean(MyBML_error_bis(n_cells,:));

JK_index_linear=mylineardecoder_error-mylineardecoder_error_bis;
JK_index_linear4=mylinear4decoder_error-mylinear4decoder_error_bis;
JK_index_popvector=mypopvector_error-mypopvector_error_bis;
JK_index_TM=myTM_error-myTM_error_bis;
JK_index_NTM=myNTM_error-MyNTM_error_bis;
JK_index_BML=myBML_error-myBML_error_bis;

end

%ici enregistrer score prop ? la configuration de population et au nbre d
%elements

end


%%
%%
%%
%%

% for i = 1:length(mypopvector_error)
% 	Mpop(i) = mean(mypopvector_error(max([1 i-10]):i)); Mtm(i) = mean(myTM_error(max([1 i-10]):i)); Mntm(i) = mean(myNTM_error(max([1 i-10]):i)); Mbml(i) = mean(myBML_error(max([1 i-10]):i));
% end
% 
% figure; hold on; plot(100-Mpop); plot(100-Mtm); plot(100.1-Mntm); plot(100-Mbml); hold off
% xlabel('Number of neurons'); ylabel('Stimulus recognition');
% title('Accuracy according to population cardinal (%)')
%  
% %ajout de 0.1 ? l'erreur du normalized TM pour visualiser deux courbes qui
% %normalement devraient etre superposees 

end

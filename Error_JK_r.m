%% ATTENTION, SCRIPT NEURONE n SPECIFIQUE, ie NEURONE n FIXE
%% D AILLEURS A TRANSFORMER EN FONCTION DE n ULTERIEUREMENT

%% SCRIPT FOR COMPUTING NEURON I DECODING IMPROVEMENT
%ici n_cells renvoie a la taille de la sous pop a laquelle l apport du neurone i est calcule
 
function [rmax,mylineardecoder_error,mylinear4decoder_error,mypopvector_error,myTM_error,myNTM_error,myBML_error,JK_index_linear,JK_index_linear4,JK_index_popvector,JK_index_TM,JK_index_NTM,JK_index_BML]=Error_JK_r(My_n, ifplots)
%My_n = neurone dont on p?se l'apport aux sous-pop dont il fait
%potentiellement partie

rmax=min(My_n-1,312-My_n);% le nb de configuration est exactement le rayon max

for r=1:rmax %rayon en nb de cellules centre (My_n) exclu

id_cells=[zeros(1,My_n-r-1),ones(1,2*r+1),zeros(1,312-My_n-r)];
id_cells_bis=id_cells;
id_cells_bis(My_n)=0;
%? changer, cette fois la composition de la population importe
%elle doit contenir le neurone i

for p=0.1:0.1:1%1 pour 10 test sets diff?rents
k=int8(p*10);
c=cvpartition(16,'Holdout',1/16);% p remplac? par 1/16 pour un one left out
id_test=test(c);
id_train=training(c);
ifplots=0;
[MyLinear_Decoder_error(r,k),MyLinear4_Decoder_error(r,k),Mypopvector_error(r,k),MyTM_error(r,k),MyNTM_error(r,k),MyBML_error(r,k)]=Decoding_DS_nc_nt(id_cells, id_train', id_test', ifplots);
[MyLinear_Decoder_error_bis(r,k),MyLinear4_Decoder_error_bis(r,k),Mypopvector_error_bis(r,k),MyTM_error_bis(r,k),MyNTM_error_bis(r,k),MyBML_error_bis(r,k)]=Decoding_DS_nc_nt(id_cells_bis, id_train', id_test', ifplots);
end
mylineardecoder_error(r)=mean(MyLinear_Decoder_error(r,:));
mylinear4decoder_error(r)=mean(MyLinear4_Decoder_error(r,:));
mypopvector_error(r)=mean(Mypopvector_error(r,:));
myTM_error(r)=mean(MyTM_error(r,:));
myNTM_error(r)=mean(MyNTM_error(r,:));
myBML_error(r)=mean(MyBML_error(r,:));

mylineardecoder_error_bis(r)=mean(MyLinear_Decoder_error_bis(r,:));
mylinear4decoder_error_bis(r)=mean(MyLinear4_Decoder_error_bis(r,:));
mypopvector_error_bis(r)=mean(Mypopvector_error_bis(r,:));
myTM_error_bis(r)=mean(MyTM_error_bis(r,:));
myNTM_error_bis(r)=mean(MyNTM_error_bis(r,:));
myBML_error_bis(r)=mean(MyBML_error_bis(r,:));

JK_index_linear=mylineardecoder_error-mylineardecoder_error_bis;
JK_index_linear4=mylinear4decoder_error-mylinear4decoder_error_bis;
JK_index_popvector=mypopvector_error-mypopvector_error_bis;
JK_index_TM=myTM_error-myTM_error_bis;
JK_index_NTM=myNTM_error-myNTM_error_bis;
JK_index_BML=myBML_error-myBML_error_bis;

end

%ici enregistrer score prop ? la configuration de population et au nbre d
%elements



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

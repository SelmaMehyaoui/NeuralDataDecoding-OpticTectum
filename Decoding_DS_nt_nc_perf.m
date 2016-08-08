clear all;
clc

%% SCRIPT FOR PLOTTING DECODERS PERFORMANCES FOR DIRECTION SELECTIVITY

MyLinear_Decoder_error=zeros(312,10);
MyLinear_Decoder_precision=zeros(312,10);
MyMatch_Decoder_error=zeros(312,10);
MyPrecision4_Decoder_error=zeros(312,10);
MyPrecision_Decoder_error=zeros(312,10);
MyMatch4_Decoder_error=zeros(312,10);
MyLinear_Decoder2_error=zeros(312,10);
MyMean4_Decoder2_error=zeros(312,10);
Mypopvector_error=zeros(312,10);
MyTM_error=zeros(312,10);
MyNTM_error=zeros(312,10);
MyBML_error=zeros(312,10);

whatdecoder=ones(8,1);

for n_cells=312:1:312 %312

for p=0.1:0.1:1 %on cross valide
k=int8(p*10);
c=cvpartition(16,'Holdout',1/16);%p
id_test=ones(1,16)';%test(c);
id_train=ones(1,16)';%training(c);
ifplots=0;

for l=1:10 %on teste 10 configurations de population diff?rentes pour ncells

id_cells=[ones(1,n_cells),zeros(1,312-n_cells)];
id_cells=id_cells(randperm(312));

saves_struct_number=0;
%saves_struct_number=l;
[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12]=Decoding_DS_nc_nt(id_cells, id_train', id_test', ifplots,whatdecoder,saves_struct_number);
MyLinear_Decoder_error(n_cells,k)=A1/10+MyLinear_Decoder_error(n_cells,k);
MyLinear_Decoder_precision(n_cells,k)=A2/10+MyLinear_Decoder_precison(n_cells,k);
MyMatch_Decoder_error(n_cells,k)=A3/10+MyMatch_Decoder_error(n_cells,k);
MyPrecision4_Decoder_error(n_cells,k)=A4/10+MyPrecision4_Decoder_error(n_cells,k);
MyPrecision_Decoder_error(n_cells,k)=A5/10+MyPrecision_Decoder_error(n_cells,k);
MyMatch4_Decoder_error(n_cells,k)=A6/10+MyMatch4_Decoder_error(n_cells,k);
MyLinear_Decoder2_error(n_cells,k)=A7/10+MyLinear_Decoder2_error(n_cells,k);
MyMean4_Decoder2_error(n_cells,k)=A8/10+MyMean4_Decoder2_error(N_cells,k);
Mypopvector_error(n_cells,k)=A9/10+Mypopvector_error(n_cells,k);
MyTM_error(n_cells,k)=A10/10+MyTM_error(n_cells,k);
MyNTM_error(n_cells,k)=A11/10+MyNTM_error(n_cells,k);
MyBML_error(n_cells,k)=A12/10+MyBML_error(n_cells,k);
end
end

end

mylineardecoder_error=mean(MyLinear_Decoder_error');
mylinear4decoder_error=mean(MyLinear4_Decoder_error');
mypopvector_error=mean(Mypopvector_error');
myTM_error=mean(MyTM_error');
myNTM_error=mean(MyNTM_error');
myBML_error=mean(MyBML_error');


for i = 1:length(mypopvector_error)
	Mlin(i)=mean(mylineardecoder_error(max([1 i-10]):i)); Mlin4(i)=mean(mylinear4decoder_error(max([1 i-10]):i)); Mpop(i) = mean(mypopvector_error(max([1 i-10]):i)); Mtm(i) = mean(myTM_error(max([1 i-10]):i)); Mntm(i) = mean(myNTM_error(max([1 i-10]):i)); Mbml(i) = mean(myBML_error(max([1 i-10]):i));
end

figure; hold on; plot(100-Mlin4); plot(100-Mpop); plot(100-Mtm); plot(100.1-Mntm); plot(100-Mbml); hold off
xlabel('Number of neurons'); ylabel('Stimulus recognition');
title('Accuracy according to population cardinal (%)')
 
%ajout de 0.1 ? l'erreur du normalized TM pour visualiser deux courbes qui
%normalement devraient etre superposees 


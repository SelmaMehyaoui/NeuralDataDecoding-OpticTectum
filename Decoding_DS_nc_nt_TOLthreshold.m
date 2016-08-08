%%% Script for fitting tolerance factor in moore-penrose pseudo-inverse
%%% computation

clc
clear all;

whatdecoder=[1 1 0 0 0 0 0 0];
id_train=[1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0];
id_test=[0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1];
ifplots=0;
saves_struct_number=0;
id_cells=ones(1,312);

for k=1:10
for l=1:10
tol1=k*0.1
tol2=l*0.1
[Linear_Decoder_Error, Linear_Decoder_Precision, Match_Decoder_Error,Precision4_Decoder_Error,Precision_Decoder_Error,Match4_Decoder_Error,]=Decoding_DS_nc_nt(id_cells, id_train, id_test, ifplots,whatdecoder,saves_struct_number,tol1,tol2);
end
end
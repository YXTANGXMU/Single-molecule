clc
clear all
close all
%Input goodtraces
[filename,filepath]=uigetfile('*.txt','Select data files','MultiSelect','off');

file=strcat(filepath,filename);
goodtraces=load (file);
 [datacut,goodtraces]=cutdata(goodtraces);

hist_draw=overlapcurves(datacut);

hist_draw1(find(hist_draw1(:,2)==-inf),2)=0;
hist_draw1(find(hist_draw1(:,2)==inf),2)=0;
hist_draw1(find(hist_draw1(:,1)==-inf),1)=0;
hist_draw1(find(hist_draw1(:,1)==inf),1)=0;

%%%one dimensional histogram
hist(hist_draw(:,2),10000);
%%%two dimensional histogram
axis_v=[-10 10 -10 10];
[MS,x,yS]=Single_Mole_Cloud(hist_draw,axis_v);





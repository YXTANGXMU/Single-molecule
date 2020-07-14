
  fs=20000; % @@@ Sampling Frequency @@@
 ft=1/fs;
 traceL=2000; %@@@@@@@@@@ traces length (how mant points per trace)@@@@@@@@@@@@@
 traceM=30000; 
lowG=-10;
 highG=1;
 lowPSD=-10;
 highPSD=1;
 n_bins=500; 
select=0; 
selectFit=0; 
cutL=500;
 prop=0.5;
 data_L=4096;

  %%cut the curves
 logG_open1=cell(1,length(logG_open));
  for i=1:length(logG_open)
       k1=min(find(logG_open{i}(1:floor(length(logG_open{i})*prop))==max(logG_open{i}(1:floor(length(logG_open{i})*prop)))));
       k2=min(find(logG_open{i}(floor(length(logG_open{i})*(1-prop)):end)==min(logG_open{i}(floor(length(logG_open{i})*(1-prop)):end))))+floor(length(logG_open{i})*(1-prop));
       logG_open1{i}=logG_open{i}(k1:min(length(logG_open{i}),k2));
 end

 high=-5;
 low=-6;
 [raw_for_flk,markA,markB]=countscut(logG_open1,high,low);

 %select the curves
Lraw=length(raw_for_flk);
z=1;
for i=1:Lraw
    if length(raw_for_flk{i}) <= traceL || length(raw_for_flk{i}) >= traceM% delte traces with lengh shorter than 500 points
        raw_for_flk{i}=[];
          z=z+1;
        continue;
    end
end
z=zeros(1,length(raw_for_flk));
 for i=1:Lraw
    z(i)=length(raw_for_flk{i});
 end
raw_for_flk=raw_for_flk(find(z~=0));
logG_open1=logG_open1(find(z~=0));

t=0;
if t== select
    Lraw=length(raw_for_flk);
    for i=1:Lraw
        raw_for_flk{i}=raw_for_flk{i}(cutL:length(raw_for_flk{i})-cutL);%%¸Ä
    end
end


% FFT and calculate PSD
for k=1:length(raw_for_flk)
    yt=transpose(raw_for_flk{k});  % Transpose the column vector to row vector
    yt=power(10,yt.*1);
    Fyt=fft(yt,data_L)/data_L*2; % FFT amplitude
    f=fs/data_L*(0:1:data_L-1); % Frequency range
    F=fs/data_L;
    A_Fyt=abs(Fyt);
    SQ_Fyt=A_Fyt.^2; % square
  
    for z=0:data_L              %  intergration strating range
        if fs/data_L*z > 100
            sp=z;
            break
        end
    end
    
    for z=sp:data_L              %  intergration end range
        if fs/data_L*z > 1000
            ep=z;
            break
        end
    end  
    int_range=(100:fs/data_L:1000);
    PSD0(k)=trapz(f(sp:ep),SQ_Fyt(sp:ep));            % @@@ Trapezoidal integration
%     @@@ correction integration  100 and 1000 are the intergation range
    if f(sp)~=100
        PSD0(k)=PSD0(k)+0.5*(((SQ_Fyt(ep+1)-SQ_Fyt(ep))/F*(1000-f(ep))+SQ_Fyt(ep))+SQ_Fyt(ep))*(1000-f(ep))-0.5*(((SQ_Fyt(sp+1)-SQ_Fyt(sp))/F*(100-f(sp))+SQ_Fyt(sp))+SQ_Fyt(sp))*(100-f(sp));
    end
    MeanG(k)=log10(mean(yt));
    PSD(k)=log10(PSD0(k)/(10^MeanG(k)));
end


meanlogG=mean(MeanG);
PSD_data= [transpose(MeanG),transpose(PSD),PSD0'];

           
[PSD_2D,X_n,Y_n]=plotPSD(PSD_data(:,1:2),highPSD,lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)

  imagesc(X_n,Y_n, PSD_2D);
set(gca,'YDir','normal')
load PSD_color
colormap(PSD_color)



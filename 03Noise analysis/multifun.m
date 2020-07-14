function funs = multifun
funs.locate=@locate;
funs.solve_equ=@solve_equ;
funs.cut_by_points=@cut_by_points;

funs.PSD_cut_1=@PSD_cut_1;
funs.draw_circle=@draw_circle;
funs.draw_could=@draw_could;
funs.PSD_cut=@PSD_cut;
funs.countsget=@countsget;
funs.defaulterrorcallback=@defaulterrorcallback;
funs.get_curve_for_tdms = @get_curve_for_tdms;
funs.cut_curves = @cut_curves;
funs.PSD_Amean = @PSD_Amean;
funs.plotPSD = @plotPSD;
funs.PSD_cut=@PSD_cut;
funs.draw_pic_1 = @draw_pic_1;
funs.PSD_matrix = @PSD_matrix;
funs.cutconductancestep=@cutconductancestep;
funs.draw_pic=@draw_pic;
funs.createFit=@createFit;
end
function   [length1]=countsget(datacut,high,low)
flag=zeros(1,size(datacut,2));
markA=zeros(1,size(datacut,2));
markB=zeros(1,size(datacut,2));
for i=1:size(datacut,2)
markB(i)=length(datacut{i});
end
length1=cell(1,size(datacut,2));
for i=1:size(datacut,2)
    for j=10:length(datacut{i})
        if datacut{i}(j)<=high&&flag(i)==0
            flag(i)=1;
           markA(i)=j;
        end
%             if datacut{i}(j,2)>high&&flag(i)==1
%             flag(i)=0;
%         end
           if datacut{i}(j)<=low&&flag(i)==1
           markB(i)=j;
           flag(i)=2;
            end
    end
end
for i=1:size(datacut,2)
    if markB(i)>markA(i)&&markA(i)~=0
    length1{i}=datacut{i}(markA(i):markB(i));
    end
end

end
function [logG_open] = cutconductancestep(data_all,additionallength1,a1,a2,b1,b2,c1,c2,d1,d2,offset,rate,V,G0,zero_set,add1)
%%%VtoI
logG_open=[];
Open=[];
data_all=data_all-offset;   %correcting voltage by offset
for i=1:size(data_all,2)
    if data_all(i)<0
        Cur0(i)=c1*data_all(i)+exp(data_all(i)*a1+b1)+d1;
    else
        Cur0(i)=c2*data_all(i)+exp(data_all(i)*a2+b2)+d2;
    end
end

for i=1:size(Cur0,2)        % Transfer I to G
    logG(i)=log10(abs(G0/((0.1/1000)/Cur0(i)-0.1)));
end


%%%截取%%%%
for i=1:size(data_all,2)-4  % Take four sampling points as
    mean_cur(i)=mean(Cur0(i:i+4));
end
k=1;
for i=1:size(data_all,2)-5
    if(mean_cur(i)<1e-7&&mean_cur(i+1)>1e-7) %find the close connecting point
        Close(k)=i;
        k=k+1;
    end
end

k=1;
for i=1:size(data_all,2)-5
    if(mean_cur(i)>1e-7&&mean_cur(i+1)<1e-7) %find the open snaping point
        Open(k)=i;
        k=k+1;
    end
end
for i=1:size(data_all,2)-20
    std_cur(i)=std(Cur0(i:i+20)); %calculate standard derivation between 20 points
end

if length(Open)~=length(Close)   % make the length of Open and Close be equal
    if length(Open)>length(Close)
        Open((length(Close)+1):end)=[];
    else
        Close((length(Open)+1):end)=[];
    end
end

if(Open(1)<Close(1)) %judge the type of the inital trace
    for i=1:length(Open)-1 %for the open initial

        Wait(i)=max(find(std_cur(Open(i):Close(i))==min(std_cur(Open(i):Close(i))))); %find the background strating point
        logG_openA{i}=logG(Open(i):(Wait(i)+Open(i)));%cut the open trace and save in the struct logG.openA
        length_open(i)=floor(Wait(i)/1000)*1000;     %why? it's just the length of Wait(i)
    end
    i=i+1;
       Wait(i)=max(find(std_cur(Open(i):end)==min(std_cur(Open(i):end)))); %find the background strating point
        logG_openA{i}=logG(Open(i):(Wait(i)+Open(i)));%cut the open trace and save in the struct logG.openA
        length_open(i)=floor(Wait(i)/1000)*1000;     %why? it's just the length of Wait(i)
    

    
else
    logG_openA=cell(1,5);
    for i=1:length(Open)-2

        Wait(i)=max(find(std_cur(Open(i):Close(i+1))==min(std_cur(Open(i):Close(i+1)))));
        logG_openA{i}=logG(Open(i):(Wait(i)+Open(i)));
        length_open(i)=floor(Wait(i)/1000)*1000;
    end
       i=i+1;
       Wait(i)=max(find(std_cur(Open(i):end)==min(std_cur(Open(i):end)))); %find the background strating point
       if (Wait(i)+Open(i))>Open(i)
        logG_openA{i}=logG(Open(i):(Wait(i)+Open(i)));%cut the open trace and save in the struct logG.openA
        length_open(i)=floor(Wait(i)/1000)*1000;     %why? it's just the length of Wait(i)
       end
    
end

%%%筛选有addtionallength
z=1;
k=1;
 logG_open=cell(1,5);
if(Open(1)<Close(1))
    
    for i=1:length(Open)-2

        if(length_open(i)+add1>additionallength1)           % extending the length of open trace
            logG_open{k}=logG(Open(i):(Wait(i)+Open(i))+(additionallength1-length_open(i)));
            k=k+1;
        end
        
        
    end
else
    for i=1:length(Open)-2

        if(length_open(i)+add1>additionallength1)
            logG_open{k}=logG(Open(i):(Wait(i)+Open(i))+(additionallength1-length_open(i)));
            k=k+1;
        end
    end
end
end


function [logG_open,logG_open_trans,test]=get_curve_for_tdms(additionallength,a1,a2,b1,b2,c1,c2,d1,d2,offset,rate,V,G0,zero_set,add1)
funs_t =TDMS_subfun;
[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else filename1{1}=filename;
end
num_file = length(filename1);
fileclass = '.tdms';
logG_open_trans=cell(1,num_file);
%  data_s=cell(1,num_file);
% test=cell(1,num_file);
parfor i = 1 : num_file
     funs_t =TDMS_subfun;
    test=funs_t.TDMS_readTDMSFile(filename1{i},filepath);
    data_s=test.data{1,3};
    funs=multifun;
    logG_open_trans{i}=funs.cutconductancestep(data_s,additionallength,a1,a2,b1,b2,c1,c2,d1,d2,offset,rate,V,G0,zero_set,add1);
    fprintf('Percentage: %03u\n',floor(i/num_file*100))
%     logG_open=[logG_open,logG_open_trans]
%     save traces_open.mat logG_open;
%     fprintf('Percentage: %03u\n',floor(i/num_file*100)); % Present the rate of progress    
end
k=zeros(1,num_file+1);
k(1)=0;
for i = 1 : num_file
     k(i+1)=k(i)+size(logG_open_trans{i},2);
end
logG_open=cell(1,k(i+1));
for i = 1 : num_file
    logG_open(1,k(i)+1:k(i+1),1)=logG_open_trans{i};
end
fprintf('Curves transform was done')
end
function [fitresult, gof] =FitPSD(X10, Y10, PSD_2D10)

[xData, yData, zData] = prepareSurfaceData( X10, Y10, PSD_2D10 );

% Set up fittype and options.
ft = fittype( 'exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
 opts.StartPoint = [0 0 0 0 0 0];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 4' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 4', 'PSD_2D10 vs. X10, Y10', 'Location', 'NorthEast' );
% % Label axes
% xlabel X10
% ylabel Y10
% zlabel PSD_2D10
% grid on
% view( -1.4, 90.0 );
end

function [Mat, x, y,G_index,bins] = plotPSD(data, max_PSD, min_PSD, max_meanG, min_meanG,n_bins)

n_bins = n_bins;
n_bins2=n_bins;
nanData = -1000;
s_C_bins = (max_PSD - min_PSD)/(n_bins);
s_D_bins = (max_meanG - min_meanG)/n_bins2;
PSD_index = linspace(max_PSD, min_PSD, n_bins);
G_index = [nanData, linspace(min_meanG, max_meanG, n_bins2)];
% Conduct_index = roundn(Conduct_index, -3);
% Dist_index = roundn(Dist_index, -3);
Matrix = zeros(n_bins, n_bins);
Matrix = [PSD_index', Matrix];
Matrix = [G_index; Matrix];
% data = roundn(data, -3);
[m, n] = size(data);

for i = 1: m
    
    x_ind = find((data(i, 1) < Matrix(1, :)+s_D_bins/2) & (data(i, 1) >= Matrix(1, :)-s_D_bins/2));
    y_ind = find((data(i, 2) < Matrix(:, 1)+s_C_bins/2) & (data(i, 2) >= Matrix(:, 1)-s_C_bins/2));
    
    Matrix(y_ind, x_ind) = Matrix(y_ind, x_ind) + 1;
    
    
end

Mat = Matrix(2:end, 2:end);
x = Matrix(1, 2:end);
y = Matrix(2:end, 1);

end


 function raw_for_flk=cut_curves(logG_open,highend,lowend)
raw_for_flk=cell(1,size(logG_open,2));
for i=1:size(logG_open,2)
    raw_for_flk{i}=logG_open{i}(find(logG_open{i}(1,:)<highend));
   raw_for_flk{i}=raw_for_flk{i}(find(raw_for_flk{i}(1,:)>lowend));
end
% save raw_for_flk.mat raw_for_flk;
% 
% figure(1)
% for i=1:10
%     
%     subplot(2,5,i);
%     n=unidrnd(length(raw_for_flk));
%     plot(raw_for_flk{n});
% end
 fprintf('Cutting curves have done')
 end
 
function  [PSD_2D_n,X_n,Y_n]=PSD_matrix(meanlogG,PSD0,MeanG,highPSD,lowPSD,highG,lowG,n_bins,PSD,n)
meanshift=(1-n)*meanlogG;
PSD_n=zeros(1,length(PSD));
for k=1:length(PSD)
    PSD_n(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD_n)];
    funs=multifun;
[PSD_2D_n,X_n,Y_n]=funs.plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
end

  function   [PSD_2D,X_n,Y_n,all_sort_1,fitresult,all_sort,No]=PSD_Amean(raw_for_flk,fs,ft,traceL,lowG,highG,lowPSD,highPSD,n_bins,select,cutL,n,flag,flagcut,Noise_upper,Noise_low,Coudu_upper,Condu_low,flagG,kkk,value1,value2)

Lraw=length(raw_for_flk);
z=1;
for i=1:Lraw
    if length(raw_for_flk{z}) <= traceL  % delte traces with lengh shorter than 500 points
        raw_for_flk(z)=[];
        continue;
    end
    z=z+1;
    if i==length(raw_for_flk)
        break;
    end
end

Lraw=length(raw_for_flk);
z=1;
for i=1:Lraw
    if length(raw_for_flk{z}) <= traceL  % delte traces with lengh shorter than 500 points
        raw_for_flk(z)=[];
        continue;
    end
    z=z+1;
    if i==length(raw_for_flk)
        break;
    end
end


t=0;
if t== select
    Lraw=length(raw_for_flk);
    for i=1:Lraw
     raw_for_flk{i}=raw_for_flk{i}(1:length(raw_for_flk{i}));
% 
% raw_for_flk{i}=raw_for_flk{i}(floor(length(raw_for_flk{i})*0.15):length(raw_for_flk{i})-floor(length(raw_for_flk{i})*0.15));
    end
end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FFT and calculate PSD
for k=1:length(raw_for_flk)
    data_L=2^nextpow2(length(raw_for_flk{k})); %   Determine Number of points (next power of 2)
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
    %int_range=(100:fs/data_L:1000);
    PSD0(k)=trapz(f(sp:ep),SQ_Fyt(sp:ep));            % @@@ Trapezoidal integration
    % @@@ correction integration  100 and 1000 are the intergation range
    if f(sp)~=100
        PSD0(k)=PSD0(k)+0.5*(((SQ_Fyt(ep+1)-SQ_Fyt(ep))/F*(1000-f(ep))+SQ_Fyt(ep))+SQ_Fyt(ep))*(1000-f(ep))-0.5*(((SQ_Fyt(sp+1)-SQ_Fyt(sp))/F*(100-f(sp))+SQ_Fyt(sp))+SQ_Fyt(sp))*(100-f(sp));
    end
    MeanG(k)=log10(mean(yt));
    PSD(k)=log10(PSD0(k)/(10^MeanG(k)));
    % PSD(k)=PSD(k)/(10^MeanG(k));
    %figure(k);
    % plot(f(1:data_L/2),SQ_Fyt(1:data_L/2));
    % xlim([90,1100]);
end
%maxpsd= max(PSD);
%minpsd= min(PSD);

PSD_2D=cell(1,length(n));
X_n=cell(1,length(n));
b_value=zeros(1,length(n));
fitresult=cell(1,length(n));
Y_n=cell(1,length(n));
meanlogG=mean(MeanG);
funs=multifun;
PSD_data= [transpose(MeanG),transpose(PSD)];
[PSD_2D{1},X_n{1},Y_n{1}]=funs.plotPSD(PSD_data,highPSD,lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)

    [fitresult{1}] =funs.createFit(X_n{1},Y_n{1}, PSD_2D{1});
    b_value(1)=fitresult{1}.b;
for i=2:length(n)
    clear PSD12
     clear PSD_data1
    meanshift=(1-n(i))*meanlogG;
for k=1:length(PSD)
    PSD12(k)=log10(PSD0(k)/(10^MeanG(k)).^n(i));
end
    PSD_data1= [transpose(MeanG),transpose(PSD12)];
[PSD_2D{i},X_n{i},Y_n{i}]=funs.plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
        [fitresult{i}] =funs.createFit(X_n{i},Y_n{i}, PSD_2D{i});
    b_value(i)=fitresult{i}.b;
end
    
% for i=2:length(n)
%     [PSD_2D{i},X_n{i},Y_n{i}]=funs.PSD_matrix(meanlogG,PSD0,MeanG,highPSD,lowPSD,highG,lowG,n_bins,PSD,n(i));
%     [fitresult{i}] =funs.createFit(X_n{i},Y_n{i}, PSD_2D{i});
%     b_value(i)=fitresult{i}.b;
% end
  No=1:1:length(n);
all_sort=[No;n;abs(b_value)]';
all_sort_1=sortrows(all_sort,3);
all_sort_2=all_sort_1(1,:);

kkk=all_sort_2(1,1);

if flagcut~=1

funs=multifun;
%  funs.draw_pic(all_sort,highG,lowG,highPSD,lowPSD,fitresult,X_n,Y_n,PSD_2D,all_sort(1,1));
 if flag==1
     clear gca
     clear gcf
     close all
     for i=1:length(No)
          funs.draw_pic(all_sort,highG,lowG,highPSD,lowPSD,fitresult,X_n,Y_n,PSD_2D,No(i));
          number=all_sort(i,2);
fileclass = '.tif';
fileclass2 = '.fig';
fileclass3='.eps';
 filename = strcat(num2str(number), fileclass);
  filename2 = strcat(num2str(number), fileclass2);
    filename3 = strcat(num2str(number), fileclass3);
 print(gcf,'-dtiff','-r600', filename)
    saveas(gcf,filename2);
        saveas(gcf,filename3);
        close all
        clear gcf
        clear gca
     end
      else
           funs.draw_pic(all_sort,highG,lowG,highPSD,lowPSD,fitresult,X_n,Y_n,PSD_2D,kkk);
          number=kkk;
fileclass = '.tif';
fileclass2 = '.fig';
fileclass3='.eps';
 filename = strcat(num2str(number), fileclass);
  filename2 = strcat(num2str(number), fileclass2);
    filename3 = strcat(num2str(number), fileclass3);
 print(gcf,'-dtiff','-r600', filename)
    saveas(gcf,filename2);
        saveas(gcf,filename3);
        close all
        clear gcf
        clear gca
     
 end
     
else
      [PSD_1, X_n_1, Y_n_1,PSD_2, X_n_2, Y_n_2]=funs.PSD_cut(PSD_2D{kkk},X_n{kkk},Y_n{kkk},Noise_upper,Noise_low,Coudu_upper,Condu_low,flagG);
    [fitresult_a] =funs.createFit(X_n_1,Y_n_1, PSD_1);
      [fitresult_b] =funs.createFit(X_n_2,Y_n_2, PSD_2);

 funs.draw_pic_1(all_sort_2,highG,lowG,highPSD,lowPSD,fitresult_a,fitresult_b,PSD_2D{kkk},X_n{kkk},Y_n{kkk},kkk,value1,value2);

hold on
end
  end
  function  draw_pic_1(all_sort,highG,lowG,highPSD,lowPSD,fitresult_a,fitresult_b,PSD_2D,X_n,Y_n,k,value1,value2)
load PSD_color.mat;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
       a = fitresult_a.a;
       b = fitresult_a.b;
       c =  fitresult_a.c;
       d = fitresult_a.d;
       e =    fitresult_a.e;
       f =   fitresult_a.f;
          a1 = fitresult_b.a;
       b1 = fitresult_b.b;
       c1 =  fitresult_b.c;
       d1 = fitresult_b.d;
       e1 =    fitresult_b.e;
       f1 =   fitresult_b.f;


%@@@@@@@@@@@@@@@@@ 拟合曲线绘制

imagesc(X_n,Y_n, PSD_2D);
set(gca,'YDir','normal')
% val_name=num2str(double(all_sort(1,2)));
%  text(0,0,val_name,'color',[R/255 G/255 B/255],'horiz','left','FontSize',45,'FontName','Arial')

colormap(PSD_color)
colorbar
hold on
     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,yRangeStart:0.01:yRangeEnd);
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f+value1);
    contour(x,y,z,5,'-k','linewidth',1.3);
    hold on
  [x1,y1]=meshgrid(xRangeStart:0.01:xRangeEnd,yRangeStart:0.01:yRangeEnd);
    z1 = exp(a1*x1.^2 + x1.*y1.*b1 + x1.*c1 + y1.*d1 + e1*y1.^2 + f1+value2);
    contour(x1,y1,z1,5,'-k','linewidth',1.3);
    hold on
zk=all_sort(k,2);
     xlabel({'$${{Conductance/G (logG_{0})}}$$'},'Interpreter','latex', 'FontSize',15')
     ylabel(['$${{Noise\, \,  Power/G (logG_{0}^n}})$$','\,\,\,','n=',num2str(zk)],'Interpreter','latex','FontSize',15')
  end
  
  function draw_pic(all_sort,highG,lowG,highPSD,lowPSD,fitresult,X_n,Y_n,PSD_2D,k)
  load PSD_color.mat;

    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
           a = fitresult{k}.a;
           b = fitresult{k}.b;
           c =  fitresult{k}.c;
           d = fitresult{k}.d;
           e =    fitresult{k}.e;
           f =   fitresult{k}.f;
       


%@@@@@@@@@@@@@@@@@ 拟合曲线绘制

imagesc(X_n{k},Y_n{k}, PSD_2D{k});
set(gca,'YDir','normal')
% val_name=num2str(double(all_sort(1,2)));
%  text(0,0,val_name,'color',[R/255 G/255 B/255],'horiz','left','FontSize',45,'FontName','Arial')

colormap(PSD_color)
colorbar
hold on
     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,yRangeStart:0.01:yRangeEnd);
     z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
    hold on
z=all_sort(k,2);
     xlabel({'$${{Conductance/G (logG_{0})}}$$'},'Interpreter','latex', 'FontSize',15')
     ylabel(['$${{Noise\, \,  Power/G (logG_{0}^n}})$$','\,\,\,','n=',num2str(z)],'Interpreter','latex','FontSize',15')
  end
  
  function [fitresult, gof] = createFit(X20, Y20, PSD_2D20)
%CREATEFIT(X20,Y20,PSD_2D20)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : X20
%      Y Input : Y20
%      Z Output: PSD_2D20
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 22-Mar-2019 15:54:41 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( X20, Y20, PSD_2D20 );

% Set up fittype and options.
ft = fittype( 'exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxFunEvals=10000;
opts.StartPoint = [0 0 0 0 0 0];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% % Create a figure for the plots.
% figure( 'Name', 'untitled fit 1' );
% 
% % Plot fit with data.
% subplot( 2, 1, 1 );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'PSD_2D20 vs. X20, Y20', 'Location', 'NorthEast' );
% % Label axes
% xlabel X20
% ylabel Y20
% zlabel PSD_2D20
% grid on
% view( 1.3, 90.0 );
% 
% % Make contour plot.
% subplot( 2, 1, 2 );
% h = plot( fitresult, [xData, yData], zData, 'Style', 'Contour' );
% legend( h, 'untitled fit 1', 'PSD_2D20 vs. X20, Y20', 'Location', 'NorthEast' );
% % Label axes
% xlabel X20
% ylabel Y20
% grid on

  end

  
  function [PSD_1, X_n_1, Y_n_1,PSD_2, X_n_2, Y_n_2]=PSD_cut(PSD_2D,X_n,Y_n,Noise_upper,Noise_low,Coudu_upper,Condu_low,flagG)
markA=1;
markB=length(X_n);
markC=1;
markD=length(Y_n);

flag=0;
for i=1:length(X_n)
        if X_n(i)>=Condu_low&&flag==0
            flag=1;
           markA=i;
        end
           if X_n(i)>=Coudu_upper&&flag==1
           markB=i;
           flag=2;
            end
end
    
flag=0;
for i=1:length(Y_n)
        if Y_n(i)<=Noise_upper&&flag==0
            flag=1;
           markC=i;
        end
           if Y_n(i)<=Noise_low&&flag==1
           markD=i;
           flag=2;
            end
end
    PSD_1=PSD_2D(markC:markD,markA:markB);
    X_n_1=X_n(markA:markB);
    Y_n_1=Y_n(markC:markD);
    if flagG==1
     PSD_2=PSD_2D(markD:length(Y_n),markA:markB);
    X_n_2=X_n(markA:markB);
    Y_n_2=Y_n(markD:length(Y_n));
    else
     PSD_2=PSD_2D(markC:markD,markB:length(X_n));
    X_n_2=X_n(markB:length(X_n));
    Y_n_2=Y_n(markC:markD);
    end
  end
  
  function   [PSD_1, X_n_1, Y_n_1,PSD_2, X_n_2, Y_n_2]=PSD_cut_1(PSD_C,X_n_C,Y_n_C,Noise_cut_upper,Noise_cut_low,Coudu_cut_upper,Condu_cut_low,flagG)
markA=1;
markB=length(X_n_C);
markC=1;
markD=length(Y_n_C);

flag=0;
for i=1:length(X_n_C)
        if X_n_C(i)>=Condu_cut_low&&flag==0
            flag=1;
           markA=i;
        end
           if X_n_C(i)>=Coudu_cut_upper&&flag==1
           markB=i;
           flag=2;
            end
end
    
flag=0;
for i=1:length(Y_n_C)
        if Y_n_C(i)<=Noise_cut_upper&&flag==0
            flag=1;
           markC=i;
        end
           if Y_n_C(i)<=Noise_cut_low&&flag==1
           markD=i;
           flag=2;
            end
end
    PSD_1=PSD_C(markC:markD,markA:markB);
    X_n_1=X_n_C(markA:markB);
    Y_n_1=Y_n_C(markC:markD);
    if flagG==1
     PSD_2=PSD_C(markD:length(Y_n_C),markA:markB);
    X_n_2=X_n_C(markA:markB);
    Y_n_2=Y_n_C(markD:length(Y_n_C));
    else
     PSD_2=PSD_C(markC:markD,markB:length(X_n_C));
    X_n_2=X_n_C(markB:length(X_n_C));
    Y_n_2=Y_n_C(markC:markD);
    end
%     for i=1:nt
%     clear a b eq1 eq2
%     syms a b eq1 eq2
%     eq1=a+b-ImB(1,i);
%     eq2=a*mt+b-ImB(mt,i);
%     [a,b]=solve(eq1,eq2,a,b);
%     aa(i)=double(a);
%     bb(i)=double(b);
  end

  
  function  draw_circle(fitresult_C,xRangeStart,xRangeEnd,yRangeStart, yRangeEnd,k1,k2)
   %%%画图
       a =fitresult_C.a;
       b =fitresult_C.b;
       c =fitresult_C.c;
       d =fitresult_C.d;
       e =fitresult_C.e;
       f =fitresult_C.f;
       %@@@@@@@@@@@@@@@@@ 拟合曲线绘制
% colorbar
% hold on
     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,yRangeStart:0.01:yRangeEnd);
    z = exp(a*(x-k1).^2 + (x-k1).*(y-k2).*b + (x-k1).*c + (y-k2).*d + e*(y-k2).^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
    hold on
  end

  
  function draw_could(all_sort,X_n,Y_n, PSD_2D,xRangeStart,xRangeEnd,yRangeStart, yRangeEnd,k)
load PSD_color
imagesc(X_n,Y_n, PSD_2D);
hold on
set(gca,'YDir','normal')
% val_name=num2str(double(all_sort(1,2)));
%  text(0,0,val_name,'color',[R/255 G/255 B/255],'horiz','left','FontSize',45,'FontName','Arial')
colormap(PSD_color)
colorbar
axis([xRangeEnd,xRangeStart,yRangeEnd, yRangeStart])
z=all_sort(k,2);
     xlabel({'$${{Conductance/G (logG_{0})}}$$'},'Interpreter','latex', 'FontSize',15')
     ylabel(['$${{Noise\, \,  Power/G (logG_{0}^n}})$$','\,\,\,','n=',num2str(z)],'Interpreter','latex','FontSize',15')
hold on
  end
  
  
  
  function    PSD_C=cut_by_points(PSD_C,points_1,points_2,points_3,points_4,X_n_C,Y_n_C)

funs = multifun;
 points_no_1=funs.locate(points_1,X_n_C,Y_n_C);
  points_no_2=funs.locate(points_2,X_n_C,Y_n_C);
   points_no_3=funs.locate(points_3,X_n_C,Y_n_C);
  points_no_4=funs.locate(points_4,X_n_C,Y_n_C);

  [k3,b3]=funs.solve_equ(points_no_1,points_no_2);  %上 
 [k4,b4]=funs.solve_equ(points_no_3,points_no_4);  %下  
 
 points_no_1=flipud(points_no_1')';
points_no_2=flipud(points_no_2')';
points_no_3=flipud(points_no_1')';
points_no_4=flipud(points_no_2')';
   [k1,b1]=funs.solve_equ(points_no_1,points_no_3);  %左   
  [k2,b2]=funs.solve_equ(points_no_2,points_no_4);  %右
 for i=1:size(PSD_C,1)
     for j=1:size(PSD_C,2)
%      if  k1*i+b1<j||k2*i+b2>j
%              PSD_C(i,j)=0;
%           end
         if  k3*j+b3>i||k4*j+b4<i
             PSD_C(i,j)=0;
          end
     end
 end
  end

  
  function points=locate(Condu,X_n_C,Y_n_C)
Coudu_value=Condu(1,1);
Noise_value=Condu(1,2);
flag_1=0;
flag_2=0;
points_1=length(Y_n_C);
points_2=length(Y_n_C);
for i=1:length(Y_n_C)
        if Y_n_C(i)<=Noise_value&&flag_1==0
            flag_1=1;
           points_2=i;
        end
           if X_n_C(i)>=Coudu_value&&flag_2==0
          flag_2=1;
           points_1=i;
            end
end
points=[points_1,points_2];
  end

  
   function [aa,bb]=solve_equ(points_1,points_3)
      syms a b eq1 eq2
      
    eq1=a*points_1(1,1)+b-points_1(1,2);
    eq2=a*points_3(1,1)+b-points_3(1,2);
    [a,b]=solve(eq1,eq2,a,b);
    
    aa=double(a);
    bb=double(b);
   if points_1(1,1)==points_3(1,1)
       aa=0;
       bb=0;
   end
    end
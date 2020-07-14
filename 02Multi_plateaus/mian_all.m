clc
clear all
close all
%Input goodtraces
[filename,filepath]=uigetfile('*.txt','Select data files','MultiSelect','off');

file=strcat(filepath,filename);
goodtraces=load (file);
 [datacut,goodtraces]=cutdata(goodtraces);

%pretreatment 
k=1;
for i=1:size(datacut,2)
    if length(datacut{i})>2000
    datacut1{k}=datacut{i};
    k=k+1;
   end
end
 k=1;
for i=1:length(datacut1)
    z=min(find(datacut1{i}(10:end,2)<-2.0));
    if isempty(z)==0
        datacut2{k}=datacut1{i};
        len_1(k)=min(find(datacut1{i}(10:end,2)<-2.0));
        datacut{k}=datacut1{i}(len_1(k)+10:max(find(datacut1{i}(:,2)>-6.0)),:);
        lengh(k)=length(datacut{k});
        k=k+1;
    end
end
datacut=datacut2;



%fragment parameters
%length 
dis=0.1;
k=10;
len_o=points(dis,datacut2{k});
kk=350000;
len=len_o;
len2=len_o;  
len1=len_o;  
%selecting condition
inte=3; 
thres1=2;    
thres2=2;    
X=1;  



%initialization

mark=zeros(floor(kk/inte),6);

for j=1:floor(kk/inte)
    mark(j,1)=inte*(j-1)+1;
     mark(j,2)=inte*(j-1)+len;
         mark(j,3)=1+mark(j,2);
     mark(j,4)=1+len1+mark(j,2);
       mark(j,5)=1+mark(j,4);
     mark(j,6)=1+len2+mark(j,4);
end

mark_c=zeros(1,length(datacut));
end_v=zeros(1,length(datacut));
meanG1=zeros(1,length(datacut));
meanG2=zeros(1,length(datacut));
mark1=cell(1,length(datacut));
mark_1=cell(1,length(datacut));

%scan the traces
l=1;
 for i=1:length(datacut)

    t1=floor(length(datacut{i})/inte)-floor(mark(1,6)/inte);
    mark1{i}=zeros(t1,3);
    mark_1{i}=zeros(t1,4);
    k=1;
 for j=1:t1
 mark1{i}(j,1)=max(datacut{i}(mark(j,1):mark(j,2),2))-min(datacut{i}(mark(j,1):mark(j,2),2));
 mark1{i}(j,2)=datacut{i}(mark(j,3),2)-datacut{i}(mark(j,4),2);
 mark1{i}(j,3)=max(datacut{i}(mark(j,5):mark(j,6),2))-min(datacut{i}(mark(j,5):mark(j,6),2));
if X==1
k1=mark1{i}(j,2)/mark1{i}(j,1);
k2=mark1{i}(j,2)/mark1{i}(j,3);
else
  k1=-mark1{i}(j,2)/mark1{i}(j,1);
k2=-mark1{i}(j,2)/mark1{i}(j,3);  
end

 if k1>thres1 && k2>thres2
     mark_1{i}(k,1)=mark(j,3);
     mark_1{i}(k,2)=k1+k2;
     mark_1{i}(k,3)=k1;
     mark_1{i}(k,4)=k2;
     k=k+1;
     label(l)=i;
          l=l+1;
 end
 end
 if k==1 
 mark_1{i}=[];
 else
   mark_1{i}=mark_1{i}(1:k-1,:);
 end
 end

%  Results
 label=unique(label);
  
 P=length(label)/length(datacut)


 function [datacut,goodtraces]=cutdata(goodtraces)

deltaConduct =0.45;
  dataSet=goodtraces;
    [m,n] = size(dataSet);
   ax=find(dataSet(:,2)==-inf);
   bx=find(dataSet(:,2)==inf);
    if(isempty(bx)==0)
   cx=[ax,bx];
    else
 cx=ax;
    end
    
   for i=3:length(cx)
      dataSet(cx(i),2)=dataSet(cx(i)-2,2);
   end
  
    %½ØÈ¡¹ì¼£
    pickIndex(1)=1;
    i=2;
    for j = 2 : m
    if abs(dataSet(j,1) - dataSet(j-1,1)) >= deltaConduct 
        pickIndex(i) = j-1;
        i=i+1;
    end
    end
    datacut=cell(1,length(pickIndex)-1);
   for i=1:length(pickIndex)-1
    datacut{i}=dataSet(pickIndex(i)+1:pickIndex(i+1),1:2);
   end
 end
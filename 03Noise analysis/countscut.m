function   [length1,markA,markB]=countscut(datacut,high,low)
flag=zeros(1,size(datacut,2));
markA=zeros(1,size(datacut,2));
markB=zeros(1,size(datacut,2));
for i=1:size(datacut,1)
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
           break
            end
    end
    if markB(i)==0
    markB(i)=length(datacut{i});
 end
end

for i=1:size(datacut,2)
    if markB(i)>markA(i)&&markA(i)~=0
    length1{i}=datacut{i}(markA(i):markB(i));
    end
end

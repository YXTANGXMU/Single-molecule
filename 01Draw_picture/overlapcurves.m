function hist_draw=overlapcurves(data_chosen)
  m=1;
for i=1:size(data_chosen,2)
m=m+size(data_chosen{i},1);
end

hist_draw=zeros(m,2);
c=1;
for i=1:length(data_chosen)
     hist_draw(c:c+size(data_chosen{i},1)-1,:)=data_chosen{i};
     c=c+size(data_chosen{i},1);
end
end

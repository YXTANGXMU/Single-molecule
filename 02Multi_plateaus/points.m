 function  p=points(dis,data)
for i=1:length(data)
if abs(abs(data(i,1)-data(2,1))-dis)<0.0001
        p=i-2;
break
end
end
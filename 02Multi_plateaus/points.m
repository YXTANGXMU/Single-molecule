 function  p=points(dis,data)
for i=1:length(data)
if abs(abs(data(i,1)-data(5,1))-dis)<0.0001
        p=i-4;
break
end
end

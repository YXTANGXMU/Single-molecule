function [markA]=curvesscan_1(datacut_single,intervals)
len_1=length(intervals);
  markA=ones(len_1+1,1);
 k=2;
 for i=1:size(datacut_single,1) 
while(k<len_1+2&&datacut_single(i,2)<=intervals(k-1))
markA(k)=i;
 k=k+1;
end
 end


end


function  [MS,x,yS]=Single_Mole_Cloud(dataoutput,axis_v)
n_bins =2500;
n_bins_1=1800;
nanData =axis_v(1);




maxConduct =axis_v(4);
 minConduct=axis_v(3);
maxDist = axis_v(2);
minDist = axis_v(1);

Conduct_index = linspace(maxConduct, minConduct, n_bins_1);
Dist_index = [ nanData linspace(minDist, maxDist, n_bins)];

kk=sortrows(dataoutput,2,'descend');
 [markA]=curvesscan_1(kk,Conduct_index);
 
  condu=cell(1,length(Conduct_index));
 for i=1:length(markA)-1
     condu{i}=kk(markA(i):markA(i+1),:);
 end
  matrix=zeros(length(Conduct_index),length(Dist_index));
  for i=1:length(condu)
 matrix(i,:)=hist(condu{i}(:,1),Dist_index);
  end
  x=Dist_index;
  y=Conduct_index;
  M= matrix;
 figure(1);
imagesc(x(2:end), -y, M);figure(gcf);
caxis([0,200])
load k.txt
 colormap(k(:,1:3));
set(gca, 'ylim', [-0.5,7.5],... 
    'YTick',(-0.5:0.5:7.5),...
    'yticklabel', {0.5;0;[]; -1.0;[]; -2.0; [];-3.0;[]; -4.0;[]; -5.0;[]; -6.0;[];-7.0;[]},'TickDir','out','xlim', [-0.5, 3.0],... 
    'XTick',(-0.5:0.5:15),'FontSize',15,'FontName','Arial')

MS=M(3:end,:);
yS=y(3:end);
end
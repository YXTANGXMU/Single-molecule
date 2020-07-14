  x=x_SMe;
  y=yS_SMe;
  M=MS_SMe;
 figure(1);
 title('2D');
imagesc(x(2:end), -y, M);figure(gcf);
caxis([0,1000])
load k.txt
 colormap(k(:,1:3));
set(gca, 'ylim', [-0.5,7.5],... 
    'YTick',(-0.5:0.5:7.5),...
    'yticklabel', {0.5;0;[]; -1.0;[]; -2.0; [];-3.0;[]; -4.0;[]; -5.0;[]; -6.0;[];-7.0;[]},'TickDir','out','xlim', [-0.5, 3.0],... 
    'XTick',(-0.5:0.5:15),'FontSize',15,'FontName','Arial')


  x=x_SAc;
  y=yS_SAc;
  M=MS_SAc;
 figure(2);
 title('2E');
imagesc(x(2:end), -y, M);figure(gcf);
caxis([0,1000])
load k.txt
 colormap(k(:,1:3));
set(gca, 'ylim', [-0.5,7.5],... 
    'YTick',(-0.5:0.5:7.5),...
    'yticklabel', {0.5;0;[]; -1.0;[]; -2.0; [];-3.0;[]; -4.0;[]; -5.0;[]; -6.0;[];-7.0;[]},'TickDir','out','xlim', [-0.5, 3.0],... 
    'XTick',(-0.5:0.5:15),'FontSize',15,'FontName','Arial')

 figure(3);
 title('3C');
  imagesc(X_n,Y_n, PSD_2D);
set(gca,'YDir','normal')
load PSD_color
colormap(PSD_color)
axis([-6.0 -4.5 -8 -5.5]) 
caxis([0,40])

  x=x_SAC1;
  y=yS_SAC1;
  M=MS_SAC1;
 figure(4);
 title('S1B');
imagesc(x(2:end), -y, M);figure(gcf);
caxis([0,1000])
load k.txt
 colormap(k(:,1:3));
set(gca, 'ylim', [-0.5,7.5],... 
    'YTick',(-0.5:0.5:7.5),...
    'yticklabel', {0.5;0;[]; -1.0;[]; -2.0; [];-3.0;[]; -4.0;[]; -5.0;[]; -6.0;[];-7.0;[]},'TickDir','out','xlim', [-0.5, 3.0],... 
    'XTick',(-0.5:0.5:15),'FontSize',15,'FontName','Arial')
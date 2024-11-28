function SURF(xdata1, ydata1, zdata1,xl,yl,zl)
colormap('jet');
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
surf(xdata1,ydata1,zdata1,'Parent',axes1,'LineStyle','none');
xlabel(xl);
ylabel(yl);
zlabel(zl);
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold','LineWidth',2);

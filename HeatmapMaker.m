function HeatmapMaker(datamat,matrixsize,colormap,lowerbound,upperbound)
xpixels = 1:matrixsize;
ypixels = 1:matrixsize;
[xmesh,ymesh]=meshgrid(xpixels,ypixels);
tb1=array2table([reshape(xmesh,[],1) reshape(ymesh,[],1) reshape(datamat,[],1)],'VariableNames',{'v1','v2','v3'});
h = heatmap(tb1,'v1','v2','ColorVariable','v3','ColorLimits',[lowerbound upperbound]);
h.Title = '';
h.XLabel = '';
h.YLabel = '';
h.ColorbarVisible = 'off';
h.Colormap = colormap;
h.ColorMethod='none';
h.MissingDataColor = [0 0 0];
h.GridVisible = 'off';
h.Position=[0 0 0.788 1];
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
Ax.CellLabelColor='none';

end
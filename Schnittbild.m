c = ShakerGeoSectioned;
m = c.createModel;
pm = PlotManager;
pm.UseFileTypeFolders = false;
pm.ExportDPI = 300;
pm.FigureSize = [1200 1000];
m.plotGeometrySetup(pm);
view(104,14);
pm.savePlots(pwd, 'Format', {'jpg','png'});
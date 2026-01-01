clc; clear;

fid = fopen('varnames.txt','r');
xname = strtrim(fgetl(fid));
yname = strtrim(fgetl(fid));
fclose(fid);

%% Load data (19x19)
corrvalue    = readmatrix('QQ_matrix.csv');
pValueMatrix = readmatrix('P_matrix.csv');

corrvalue_flip_r=flip(corrvalue,1);
corrvalue_flip_c=flip(corrvalue_flip_r,2);

pValueMatrix_flip_r=flip(pValueMatrix,1);
pValueMatrix_flip_c=flip(pValueMatrix_flip_r,2);

corrvalue    = corrvalue_flip_c(1:19,1:19);
pValueMatrix = pValueMatrix_flip_c(1:19,1:19);

pValueThreshold1 = 0.01;
pValueThreshold2 = 0.05;
pValueThreshold3 = 0.1;
imagesc(corrvalue);
colormap(summer);
colorbar;
[numRows, numCols] = size(corrvalue);
for row = 1:numRows
for col = numCols:-1:1
if pValueMatrix(row, col) < pValueThreshold1
pValueStr = '***';
text(col, row, pValueStr, 'HorizontalAlignment', 'center');
elseif pValueMatrix(row, col) > pValueThreshold3
pValueStr = '';
text(col, row, pValueStr, 'HorizontalAlignment', 'center');
elseif pValueThreshold2 < pValueMatrix(row, col) < pValueThreshold3
pValueStr = '**';
text(col, row, pValueStr, 'HorizontalAlignment', 'center');
else
pValueStr = '*';
text(col, row, pValueStr, 'HorizontalAlignment', 'center');
end
end
end
title(['Cross-Quantile Correlation Heatmap between ' yname ' and ' xname]);
xlabel(['Quantiles of ' xname]);
ylabel(['Quantiles of ' yname]);
ax = gca;
ax.XAxis.TickValues = 1:numRows;
ax.YAxis.TickValues = 1:numCols;
xTickLabels = ax.XAxis.TickLabels;
for i = 1:length(xTickLabels)
xTickLabels{i} = sprintf('%0.2f',i * 0.05);
end
ax.XAxis.TickLabels = xTickLabels;
yTickLabels = ax.YAxis.TickLabels;
for i = 1:length(yTickLabels)
yTickLabels{i} = sprintf('%0.2f',(length(yTickLabels) - i + 1) * 0.05);
end
ax.YAxis.TickLabels = yTickLabels;
set(gcf,'WindowState', 'Fullscreen');
saveas(gcf,'Cross-Quantile Correlation Heatmap.png');

if exist('varnames.txt','file')
    delete('varnames.txt');
end
if exist('QQ_matrix.csv','file')
    delete('QQ_matrix.csv');
end
if exist('P_matrix.csv','file')
    delete('P_matrix.csv');
end
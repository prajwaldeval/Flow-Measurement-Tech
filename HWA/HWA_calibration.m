clc; clear variables; close all; 

%%  Config
urange = 0:2:20; 
FileRoot = 'C:\Users\LocalAdmin\Documents\FMT\2022\G20\Calibration_'; 

for ii = 1:numel(urange)
%     fileIn = sprintf(fileName,urange(ii))
if urange(ii) < 10
    fileIn = [FileRoot '00' num2str(urange(ii),'%3d')];
else
        fileIn = [FileRoot '0' num2str(urange(ii),'%3d')];
end
    delimiter = ' ';
    startRow = 23;
    formatSpec = '%s';
    try
        fileID = fopen(fileIn,'r');
    catch
        fileID = fopen(fileIn{1},'r');
    end
    tmp = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    tmp2 = strrep(tmp{1},',','.');
    
    fileOut = strcat(fileIn,'_mat');
    fileId = fopen(fileOut,'w');
    fprintf(fileId, '%s\n', tmp2);
    fclose(fileID);
    
    
    data  = dlmread(fileOut); 
    u(ii) = mean(data(:,2));
%     pause 
end

%%  Fit Poly
x = u; 
xplot = linspace(x(1),x(end),100); 
y = urange; 
c = polyfit(x,y,4); 

figure(1); clf; set(gcf,'color','w','position',[381 558 859 420]); 
plot(y,x,'k.','markersize',10); 
hold on; 
plot(polyval(c,xplot),xplot)
xlabel('u [m/s]'); ylabel('E [V]'); 
set(gca,'fontsize',14)

% dim = [.2 .5 .3 .3];
% str = 'u = (%2.2f) E^4 + (%2.2f) E^3 + (%2.2f) E^2 + (%2.2f) E + (%2.2f)';
% annotation('textbox',dim,'String',sprintf(str,c),'FitBoxToText','on');
%title(sprintf(str,c))

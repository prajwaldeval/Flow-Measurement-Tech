clc
clear
close all

%% Data Extraction 

AoA = '_00'; %%%%%%%%%%%%%%%%  CHOOSE ANGLE OF ATTACK  %%%%%%%%%%%%%%%%%%%%

% Calibration data
FileRoot = 'C:\Users\prajw\Documents\Study Stuff\MSc-1\FlowMeasurementTech\FMT\HWA\G20\Calibration_'; 


time = [];
voltage = [];
samples = [];
sampling_rate = [];
mean_voltage = [];
sigma_voltage = [];
data_cali = [];  % Calibration data

for i = 0:2:20
    
if i < 10
    fileIn = [FileRoot, '00', num2str(i,'%3d')];
else
    fileIn = [FileRoot, '0', num2str(i,'%3d')];
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

data_1 = dlmread(fileOut);

data_cali = [data_cali data_1];

% time = [time  data_1( : , 1 )] ;
% voltage = [voltage  data_1( : , 2 )] ;
% samples = [samples  numel( voltage(:,i/2 +1 ))] ;
% sampling_rate = [sampling_rate  samples(:,i/2 +1 )./max( time )] ;
% mean_voltage = [mean_voltage  mean( voltage(:,i/2 +1 ) )] ;
end

for i = 1:11
time = [time data_cali( : , 2*i -1 )] ;
voltage = [voltage data_cali( : , 2*i )] ;
mean_voltage = [mean_voltage mean(voltage(:,i))] ;
sigma_voltage = [sigma_voltage std(voltage(:,i))] ;
end

% Measured data
FileRoot = 'C:\Users\prajw\Documents\Study Stuff\MSc-1\FlowMeasurementTech\FMT\HWA\G20\Measurement_';

data_meas = [];  % Measurement data
time_meas = [];
voltage_meas = [];
mean_voltage_meas = [];
sigma_voltage_meas = [];
samples_t = [];
sampling_rate_t = [];

for i = 30:4:106

    fileIn = [FileRoot num2str(i,'%3d') AoA];
     
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

    data_meas = [data_meas dlmread(fileOut)];
end

for i = 1:20
time_meas = [time_meas data_meas(:,2*i-1)] ;
voltage_meas = [voltage_meas data_meas( : , 2*i )] ;
mean_voltage_meas = [mean_voltage_meas mean(voltage_meas(:,i))] ;
sigma_voltage_meas = [sigma_voltage_meas std(voltage_meas(:,i))] ;
samples_t = [samples_t length(voltage_meas( : , i))] ;
sampling_rate_t = [sampling_rate_t samples_t( : , i)./(time_meas( end , i) - time_meas(1 , i) )] ;
end

%% Outliers removal outside +-3σ
% Calibration data outliers
for i = 1:11
    for j = 1:10000
        if abs(voltage(j,i)-mean_voltage(i))>3*sigma_voltage(i)
            voltage(j,i) = mean_voltage(i);
        end
    end
    mean_voltage(i) = mean(voltage(:,i)) ;
    sigma_voltage(i) = std(voltage(:,i)) ;
end

% Measured data outliers 
for i = 1:20
    for j = 1:850
        if abs(voltage_meas(j,i)-mean_voltage_meas(i))>3*sigma_voltage_meas(i)
            voltage_meas(j,i) = mean_voltage_meas(i);
        end
    end
    mean_voltage_meas(i) = mean(voltage_meas(:,i)) ;
    sigma_voltage_meas(i) = std(voltage_meas(:,i)) ;
end

%% Calibration 
velocity = 0:2:20 ;
kings_coefficients = polyfit(mean_voltage , velocity , 4) ;
voltage_range = linspace(mean_voltage (1) , mean_voltage (end) ) ;
p = polyval( kings_coefficients , voltage_range ) ;

% figure (1)
% plot ( velocity, mean_voltage , ' o ' )
% hold on
% plot (p , voltage_range ) ;
% title( ' Calibration Curve ' )
% xlabel ( 'U [m/s ] ' )
% ylabel ( 'E [V] ' )
% xlim ( [ 0 20])
% grid
% hold off

%% Velocity 
c = kings_coefficients;
c1 = c(5) ;
c2 = c(4) ;
c3 = c(3) ;
c4 = c(2) ;
c5 = c(1) ;

velocity_meas = c5*( voltage_meas.^4) + c4*( voltage_meas.^3) + c3*( voltage_meas.^2) + c2*voltage_meas + c1;
mean_velocity_meas = mean( velocity_meas ) ;
sigma_velocity_meas = std( velocity_meas ) ;

%% Time Series
N = 3; % Number of measurement for autocorrelation

% figure (2)
% scatter(time(:,N) , voltage(:,N), 5,'filled' )
% hold on
% plot (time(:,N) , (mean_voltage(N)*ones(1 , length(time(:,N)) ) ) , ' r ' )
% hold on
% plot (time(:,N) , (mean_voltage(N) + 3*sigma_voltage(N) )*ones(1 , length(time(:,N)) ) , 'k ' )
% hold on
% plot (time(:,N) , (mean_voltage(N) - 3*sigma_voltage(N) )*ones(1 , length(time(:,N)) ) , 'k ' )
% title( 'Time Series of Signal ' )
% xlabel( 'T [s] ' )
% ylabel( 'E [V] ' )
% legend( 'E' , 'E_{mean} ' , '3 \sigma_{E} ' )
% grid
% hold off

% figure (3)
% plot ((1000*time(1:300,N) ) , voltage(1:300,N) )
% hold on
% plot ((1000*time(1:300,N) ) , mean_voltage(N)*ones (1 , length (time(1:300,N) ) ) , ' r ' )
% hold on
% plot ((1000*time(1:300,N) ) , (mean_voltage(N) + sigma_voltage(N))*ones(1,length(time(1:300,N))),'k')
% hold on
% plot ((1000*time(1:300,N) ) , (mean_voltage(N) - sigma_voltage(N))*ones(1,length(time(1:300,N))),'k')
% title( 'Time Series of Signal ' )
% xlabel ( 'T [ms] ' )
% ylabel ( 'E [V] ' )
% % ylim ( [ 1.7  1.74 ] )
% legend ( 'E' , 'E_{mean} ' , ' \sigma_{E} ' )
% grid
% hold off

%% Correlation data
FileRoot = 'C:\Users\prajw\Documents\Study Stuff\MSc-1\FlowMeasurementTech\FMT\HWA\G20\CorrelationTest'; 

time_cor = [];
voltage_cor = [];
samples_cor = [];
sampling_rate_cor = [];
mean_voltage_cor = [];
sigma_voltage_cor = [];
data_cor = [];  % Correlation data


fileIn = FileRoot;

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

data_cor = dlmread(fileOut);


% time = [time  data_1( : , 1 )] ;
% voltage = [voltage  data_1( : , 2 )] ;
% samples = [samples  numel( voltage(:,i/2 +1 ))] ;
% sampling_rate = [sampling_rate  samples(:,i/2 +1 )./max( time )] ;
% mean_voltage = [mean_voltage  mean( voltage(:,i/2 +1 ) )] ;



time_cor = data_cor(:,1) ;
voltage_cor = data_cor(:,2) ;
mean_voltage_cor = mean(voltage_cor) ;
sigma_voltage_cor = std(voltage_cor) ;


%% Autocorrelation
conv_crit = 0.1; %%%%%%%%%%%%%%%%%% Choose %%%%%%%%%%%%%%%%%%%%%
dt = time(2) - time(1) ;

numero = 1000;

deviations_volt = voltage(:,N) - mean_voltage(N)  ;
acorr_volt = xcorr( deviations_volt ) ; % , 'coeff' , 'coeff'
deviations_vel = velocity_meas(:,N) - mean_velocity_meas(N);
acorr_vel = xcorr( deviations_vel , 'coeff') ; 
deviations_cor = voltage_cor - mean_voltage_cor;
acorr_cor = xcorr( deviations_cor , 'coeff') ; 

acorr = acorr_cor;
rho = acorr( length ( voltage_cor ) : ( length ( voltage_cor ) + numero) ) ;

T = ( 0 : ( length ( rho ) -1) ) .* dt ;
freq = 1./(T) ;
% figure (4)
% plot ((T*1000) , rho , 'k ' )
% title( ' Autocorrelation ' )
% xlabel ( 'T [ms] ' )
% ylabel ( ' \rho [ - ] ' )
% grid
% hold on
% plot ((T*1000) , conv_crit*ones (1 , length (T) ) , ' r ' )
% legend ( 'Correlation Coefficient ' , ' Convergence Criteria ' )
% hold off
% figure (5)
% semilogx ( freq , rho , 'k ' )
% title( ' Autocorrelation ' )
% xlabel ( ' f [Hz ] ' )
% ylabel ( ' \rho [ - ] ' )
% grid
% hold on
% plot ( freq , conv_crit*ones (1 , length (T) ) , ' r ' )
% legend ( ' Correlation Coefficient ' , ' Convergence Criterion ' , ' location ' , ' northwest ' )
% hold off

%%%%%%%%%%%%%%%%%%%%%%
% nt = 10000;
% r = xcorr(deviations); r = r(nt:end); r = r/r(1);
% figure(10); clf; set(gcf,'color','w','position',[662 556 560 420]); 
% subplot(211); 
% plot(time(1:nt,N)*1000,r)
% xlabel('time [ms]'); ylabel('\rho');
% set(gca,'fontsize',14,'xlim',[0 50]); grid on; box on; 
% % highlight first zero crossing
% indz = find(r <= 0); indz = indz(1); 
% dtplot = abs(r(indz))/abs((r(indz-1)-r(indz))/dt);
% fi = 1/(time(indz,N)-dtplot); 
% % if (abs(r(indz)) > abs(r(indz-1)))
% %     indz = indz-1; 
% % end
% hold on; plot((time(indz,N)-dtplot)*1000,0,'.','markersize',20);
% txtstr = sprintf('$$[T_I = %1.2f ms]$$',(time(indz,N)-dtplot)*1000);
% ht = annotation('textbox',[.15 .85 0.4 0.03],'String',txtstr,'FitBoxToText','off',...
%                 'edgecolor','none','fontsize',15,'HorizontalAlignment','left','color','k','interpreter','latex');
%             
% subplot(212); 
% semilogx(1./time(1:nt,N),r)
% xlabel('frequency [Hz]'); ylabel('\rho'); 
% set(gca,'fontsize',14,'xlim',[1e1 1e4]); grid on; box on; 
% hold on; semilogx(1/(time(indz,N)-dtplot),0,'.','markersize',20);
% txtstr = sprintf('[f = %4i Hz]',round(1/(time(indz,N)-dtplot)));
% ht = annotation('textbox',[.69 .2 0.4 0.03],'String',txtstr,'FitBoxToText','off',...
%                 'edgecolor','none','fontsize',15,'HorizontalAlignment','left','color','k','interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%

%% Frequency Domain Analysis - Spectral Analysis
na = 7; % Averaging coefficient
ov = 0; % Number of samples for overlap
SR = 10000; % Sampling frequency [Hz]
w = hann( floor( length ( voltage(:,3) )/na) ) ; % Window
w_ov = length(w)*ov ; % Window overlap
[ psd , f ] = pwelch( voltage(:,3), w, w_ov, [ ] , SR) ;
f_bin = f (2) - f (1) ; % Frequency bin width [Hz ]
PSD = psd*f_bin; % Power spectral density [ unit /Hz ]
% figure(6)
% loglog( f , PSD)
% title( ' Spectral Analysis ' )
% xlabel( ' f [Hz ] ' )
% ylabel( 'PSD [m^2/ s^2/Hz ] ' )
% xlim([10 5000])
% grid

%% Amplitude Analysis


% Time Series of Velocity
MN = 20; % Measurement number
% figure(7)
% P(1) = plot( time_meas( : ,MN) , velocity_meas( : ,MN) ) ;
% hold on
% P(2) = plot( time_meas( : ,MN) ,( mean_velocity_meas(MN)*ones(1 , length( time_meas( : ,MN)))), ' r ' ) ;
% hold on
% P(3) = plot(time_meas(:,MN),(mean_velocity_meas(MN)+sigma_velocity_meas(MN))*ones(1,length(time_meas(:,MN))),'k');
% hold on
% P(4) = plot( time_meas( : ,MN) ,( mean_velocity_meas(MN) - sigma_velocity_meas(MN) )*ones(1 , length( time_meas( : ,MN) ) ) , 'k ' ) ;
% hold on
% P(5) = plot( time_meas( : ,MN) ,( mean_velocity_meas(MN)+3*sigma_velocity_meas(MN) )*ones(1 , length( time_meas ( : ,MN) ) ) , ' --k ' ) ;
% hold on
% plot(time_meas(:,MN),(mean_velocity_meas(MN)-3*sigma_velocity_meas(MN))*ones(1,length(time_meas(:,MN))),'--k');
% title( 'Time Series of Velocity ' )
% xlabel( 'T [ s ] ' )
% ylabel( 'U [m/s ] ' )
% ylim( [ 0 20])
% legend(P( [ 1 2 3 5 ] ) , 'U' , 'U_{mean} ' , 'U_{RMS} ' , ' 3\sigma_{U} ' , ' location ' , ' northwest ' )
% grid
% hold off
%% Velocity Profile
Y = 34:4:106;
% figure (8)
% plot( mean_velocity_meas(2:20) , Y, ' -x r ' )
% % hold on
% % plot( mean_velocity_meas(21:40) , Y, ' -x b ' )
% % hold on
% % plot( mean_velocity_meas(41:60) , Y, ' -x k ' )
% title( ' Velocity Profile' )
% xlabel( 'U_{mean} [m/s ] ' )
% ylabel( 'Y [mm] ' )
% grid
% % legend ( ' \alpha = 0\circ ' , ' \alpha = 5\circ ' , ' \alpha = 15\circ ' , ' location ' , ' east ' )
% hold off
% figure(9)
% plot( sigma_velocity_meas(2:20) , Y, ' -x r ' )
% % hold on
% % plot( sigma_velocity_meas(22:42) , Y, ' -x b ' )
% % hold on
% % plot( sigma_velocity_meas(43:63) , Y, ' -x k ' )
% title( ' Velocity Profile' )
% xlabel( 'U_{RMS} [m/s]' )
% ylabel( 'Y [mm]' )
% grid
% % legend ('\alpha = 0\circ', '\alpha = 5\circ', '\alpha = 15\circ')
% hold off

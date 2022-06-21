clc; clear; close all;



%%
option = 3; 
if (option == 3)
    fileIn = 'C:\Users\mvret\Desktop\Flow-Measurement-Tech\HWA\G20\Calibration_004';
    data  = dlmread(fileIn,'',23,0);
    t     = data(:,1); 
    E     = data(:,2); 
    u = 46.256*E.^4 - 226.55*E.^3 + 444.82*E.^2 - 400.82*E + 136.52 ;
    clear data fileIn
    yrange = [-.2 .2];
    ytickvar  = [-.1:.1:.1];
    ystr = 'E'' [V]';
    nsub = 850;
end
dt = t(2)-t(1); 
config = 2;
nt = numel(t);
%% time series 
figure(1); clf; set(gcf,'color','w','position',[70 556 560 420]); 
subplot(211)
plot(t(1:nt),u(1:nt)-mean(u(1:nt)));
grid on;  xlabel('t [s]'); box on; ylabel(ystr);
set(gca,'fontsize',14,'ylim',yrange,'ytick',ytickvar); 
hold on; plot(t(nsub)*1000*[1 1],[floor(min(u)) ceil(max(u))])

subplot(212)
plot(t(nsub)*1000,u(nsub)-mean(u(1:nsub)));
grid on;  xlabel('t [ms]'); box on; ylabel(ystr);
set(gca,'fontsize',14,'ylim',yrange,'ytick',ytickvar);  
hold on; 
%% statistics
if (config == 1)
    u_sub  = u; 
    u_mean = mean(u);
    u_std = std(u);
    up = u - u_mean;
elseif (config == 2)
    u_sub  = u(1:nt);
    u_mean = mean(u(1:nt));
    u_std  = std(u(1:nt));
    up     = u - mean(u);%u_mean;
    upp    = u_sub - u_mean;
end
% Add mean to plot
subplot(211); hold on; 
% plot([t(1) t(nt)],u_mean*[1 1],'k--','linewidth',1);
set(gca,'fontsize',14); 

%% Autocorrelation
r = xcorr(upp); r = r(nt:end); r = r/r(1);
figure(2); clf; set(gcf,'color','w','position',[662 556 560 420]); 
subplot(211); 
plot(t(1:nt)*1000,r)
xlabel('time [ms]'); ylabel('\rho');
set(gca,'fontsize',14,'xlim',[0 50]); grid on; box on; 
% highlight first zero crossing
indz = find(r <= 0); indz = indz(1); 
dtplot = abs(r(indz))/abs((r(indz-1)-r(indz))/dt);
fi = 1/(t(indz)-dtplot); 
% if (abs(r(indz)) > abs(r(indz-1)))
%     indz = indz-1; 
% end
hold on; plot((t(indz)-dtplot)*1000,0,'.','markersize',20);
txtstr = sprintf('$$[T_I = %1.2f ms]$$',(t(indz)-dtplot)*1000);
ht = annotation('textbox',[.15 .85 0.4 0.03],'String',txtstr,'FitBoxToText','off',...
                'edgecolor','none','fontsize',15,'HorizontalAlignment','left','color','k','interpreter','latex');
            
subplot(212); 
semilogx(1./t(1:nt),r)
xlabel('frequency [Hz]'); ylabel('\rho'); 
set(gca,'fontsize',14,'xlim',[1e1 1e4]); grid on; box on; 
hold on; semilogx(1/(t(indz)-dtplot),0,'.','markersize',20);
txtstr = sprintf('[f = %4i Hz]',round(1/(t(indz)-dtplot)));
ht = annotation('textbox',[.69 .2 0.4 0.03],'String',txtstr,'FitBoxToText','off',...
                'edgecolor','none','fontsize',15,'HorizontalAlignment','left','color','k','interpreter','latex');
% return
%% power spectrum
[Pxx,F] = pwelch(X,WINDOW,NOVERLAP,NFFT,Fs)
f_acq = (t(2)-t(1))^-1;
L = 10000;
win = hann(L);
Noverlap = L/2;
[Pxx,freq] = pwelch(up,win,Noverlap,L,f_acq);
figure(3); clf; set(gcf,'color','w','position',[1253 556 560 420]); 
loglog(freq,Pxx), grid on, xlabel('f [Hz]'), ylabel('PSD [m^2/s^2/Hz]')

indz2 = find(Pxx == max(Pxx)); 
hold on; loglog(freq(indz2),Pxx(indz2),'.','markersize',20);
txtstr = sprintf('[f = %4i Hz]',round(freq(indz2)));
ht = annotation('textbox',[.61 .895 0.4 0.03],'String',txtstr,'FitBoxToText','off',...
                'edgecolor','none','fontsize',15,'HorizontalAlignment','left','color','k','interpreter','latex');
set(gca,'fontsize',14,'xlim',[10 5e3])

function u = kingslaw(E)
    kingslaw = 46.256*E.^4 - 226.55*E.^3 + 444.82*E.^2 - 400.82*E + 136.52;
end
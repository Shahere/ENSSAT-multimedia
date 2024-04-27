% Load the IRC file
load IRC_1002_C_HRIR.mat
l_eq_hrir_S
%%
Fs = l_eq_hrir_S.sampling_hz;
N = 73;
l_eq_hrir_S.azim_v(N);
l_eq_hrir_S.elev_v(N);
RepImpLeft = l_eq_hrir_S.content_m(N,:);
RepImpRight = r_eq_hrir_S.content_m(N,:);
HLeft = fft(RepImpLeft);
HRight = fft(RepImpRight);
Indices = (1:1+length(RepImpRight)/2);
F = (Indices - 1) * Fs / length(RepImpRight);
figure;
semilogx(F, 20*log10(abs(HLeft(Indices))));
hold on;
semilogx(F, 20*log10(abs(HRight(Indices))), 'm');
grid;
axis([100 Fs/2 -30 10]);
xlabel('Frequence (Hz)');
ylabel('Amplitude (dB)');
%%
N = 79;
l_eq_hrir_S.azim_v(N);
l_eq_hrir_S.elev_v(N);
RepImpLeft = l_eq_hrir_S.content_m(N,:);
RepImpRight = r_eq_hrir_S.content_m(N,:);
HLeft = fft(RepImpLeft);
HRight = fft(RepImpRight);
Indices = (1:1+length(RepImpRight)/2);
F = (Indices - 1) * Fs / length(RepImpRight);
figure;
semilogx(F, 20*log10(abs(HLeft(Indices))));
hold on;
semilogx(F, 20*log10(abs(HRight(Indices))), 'm');
grid;
axis([100 Fs/2 -30 10]);
xlabel('Frequence (Hz)');
ylabel('Amplitude (dB)');
%% Angle azimutal de 0° et angles d'élévations
figure;
for n=0:6
    N = 1 + (24 * n);
    l_eq_hrir_S.azim_v(N);
    l_eq_hrir_S.elev_v(N);
    RepImpLeft = l_eq_hrir_S.content_m(N,:);
    RepImpRight = r_eq_hrir_S.content_m(N,:);
    HLeft = fft(RepImpLeft);
    HRight = fft(RepImpRight);
    Indices = (1:1+length(RepImpRight)/2);
    F = (Indices - 1) * Fs / length(RepImpRight);
    semilogx(F, 20*log10(abs(HLeft(Indices))));
    hold on;
    semilogx(F, 20*log10(abs(HRight(Indices))), 'm');
end
legend('L1','R1','L25','R25','L49','R49','L73','R73','L97','R97','L121','R121','L145','R145');
grid;
axis([100 Fs/2 -30 10]);
xlabel('Frequence (Hz)');
ylabel('Amplitude (dB)');
%% Angle azimutal de 90° et angles d'élévations
figure;
for n=0:6
    N = 7 + (24 * n);
    l_eq_hrir_S.azim_v(N);
    l_eq_hrir_S.elev_v(N);
    RepImpLeft = l_eq_hrir_S.content_m(N,:);
    RepImpRight = r_eq_hrir_S.content_m(N,:);
    HLeft = fft(RepImpLeft);
    HRight = fft(RepImpRight);
    Indices = (1:1+length(RepImpRight)/2);
    F = (Indices - 1) * Fs / length(RepImpRight);
    semilogx(F, 20*log10(abs(HLeft(Indices))));
    hold on;
    semilogx(F, 20*log10(abs(HRight(Indices))), 'm');
end
legend('L1','R1','L25','R25','L49','R49','L73','R73','L97','R97','L121','R121','L145','R145');
grid;
axis([100 Fs/2 -30 10]);
xlabel('Frequence (Hz)');
ylabel('Amplitude (dB)');

%% HRTF
Fs = 44100;
duree = 4;
t=transpose(linspace(0, duree, duree*Fs));
bruit = rand(size(t));
Sig44 = floor(bruit+0.001);
%soundsc(Sig44, Fs)

%Filtrage du signal audio par 2 reponses impulsionnelles

Zin = [];
[Zout, OutLeft] = filtrage(RepImpLeft, 1, Sig44, Zin);
Zin = [];
[Zout, OutRight] = filtrage(RepImpRight, 1, Sig44, Zin);
Out(:,1) = OutLeft;
Out(:,2) = OutRight;
audiowrite('SonSpatialiser.wav', Out, Fs)
soundsc(Out, Fs)
%wavwrite(Out, Fs, 'SonSpatialiser.wav');

%% barber
[data, Fs] = audioread('BarberShop.mp3');
soundsc(data, Fs);

%% ITD ILD
p.Fs = 44100;
p.dur = 4;
t=linspace(0, p.dur, p.dur*p.Fs);
centerFreq=440;
width=200;
y=randn(size(t));
F=complex2real(fft(y), t);
Gaussian = exp(-(F.freq-centerFreq).^2/width^2)';
plot(F.freq, Gaussian)
F.amp = F.amp.*transpose(Gaussian);
yy = transpose(real(ifft(real2complex(F))));
sound(yy, p.Fs);

%% Graphique de x et y en fonction du temps
figure;
p.Fs = 44100;
p.dur = 4;
t=linspace(0, p.dur, p.dur*p.Fs);
x = [];
y = [];
dl = [];
dr = [];
for n=1:length(t)
    x(n) = -50 + (100 * n / length(t));
    y(n) = 10;
    dl(n) = sqrt((x(n) - (17.5 / 2))^2 + 10^2);
    dr(n) = sqrt((x(n) + (17.5 / 2))^2 + 10^2);
end
plot(t, x, t, y, t, dl, t, dr);
legend('x', 'y', 'dLeft', 'dRight');
xlabel('temps');
ylabel('Distance');

%% temps de trajet
figure;
c = 345;
p.Fs = 44100;
p.dur = 4;
t=linspace(0, p.dur, p.dur*p.Fs);
x = [];
y = [];
dl = [];
dr = [];
tl = [];
tr = [];
for n=1:length(t)
    x(n) = -50 + (100 * n / length(t));
    y(n) = 10;
    dl(n) = sqrt((x(n) - (17.5 / 2))^2 + 10^2);
    dr(n) = sqrt((x(n) + (17.5 / 2))^2 + 10^2);
    tl(n) = dl(n) / c;
    tr(n) = dr(n) / c;
end
plot(t, tl, t, tr);
legend('tLeft', 'tRight');
xlabel('Temps de trajet');
ylabel('Temps');

%% Effet doppler
figure;
p.Fs = 44100;
p.dur = 4;
t=linspace(0, p.dur, p.dur*p.Fs);
centerFreq=440;
width=200;
yt=randn(size(t));
F=complex2real(fft(yt), t);
Gaussian = exp(-(F.freq-centerFreq).^2/width^2)';
F.amp = F.amp.*transpose(Gaussian);
yy = transpose(real(ifft(real2complex(F))));

c = 345;
x = [];
y = [];
dl = [];
dr = [];
tl = [];
tr = [];
for n=1:length(t)
    x(n) = -50 + (100 * n / length(t));
    y(n) = 10;
    dl(n) = sqrt((x(n) - (17.5 / 2))^2 + 10^2);
    dr(n) = sqrt((x(n) + (17.5 / 2))^2 + 10^2);
    tl(n) = dl(n) / c;
    tr(n) = dr(n) / c;
end
sLevelL = yy./(dl.^2)';
sLevelR = yy./(dr.^2)';
yy = [sLevelL, sLevelR];
yy = yy/max(yy(:));
soundsc(yy, p.Fs)

dt = 1/(t(2)-t(1));
speedL = diff(sqrt(dl))*dt;
speedR = diff(sqrt(dr))*dt;
time = t(1:end-1);
plot(time, speedL, time, speedR)
legend('speedL', 'speedR');
xlabel('Temps');
ylabel('Vitesse');

%% Stereo
p.Fs = 44100;
p.dur = 4;
t=linspace(0, p.dur, p.dur*p.Fs);
centerFreq=440;
width=200;
yt=randn(size(t));
F=complex2real(fft(yt), t);
Gaussian = exp(-(F.freq-centerFreq).^2/width^2)';
F.amp = F.amp.*transpose(Gaussian);
yy = transpose(real(ifft(real2complex(F))));

c = 345;
x = [];
y = [];
dl = [];
dr = [];
tl = [];
tr = [];
for n=1:length(t)
    x(n) = -50 + (100 * n / length(t));
    y(n) = 10;
    dl(n) = sqrt((x(n) - (17.5 / 2))^2 + 10^2);
    dr(n) = sqrt((x(n) + (17.5 / 2))^2 + 10^2);
    tl(n) = dl(n) / c;
    tr(n) = dr(n) / c; % tauR et tauL
end
%sound level
sLevelL = yy./(dl.^2)';
sLevelR = yy./(dr.^2)';
yy = [sLevelL, sLevelR];
yy = yy/max(yy(:));
%effet doppler
dt = 1/(t(2)-t(1));
speedL = diff(sqrt(dl))*dt;
speedR = diff(sqrt(dr))*dt;
%%spacialisé
yy = zeros(length(t), 2);
yy(:,2) = interp1(t, sLevelL, t-tl);
yy(:,1) = interp1(t, sLevelR, t-tr);
yy(isnan(yy)) = 0;
yy = yy/max(yy(:));
soundsc(yy,p.Fs);

plot(t, yy(:,1), t, yy(:,2));
legend('Left', 'Right');
xlabel('Temps');
ylabel('Amplitude');

%% espace 3D
load handel
s = y;
p.Fs = Fs;
p.dur = length(s)/p.Fs;
t=transpose(linspace(0,p.dur, length(s)));
x = 20*sin(2*pi*t/p.dur);
y = 5*cos(2*pi*t/p.dur);

[yy, p]=animateAuditoryMotion(p,s,x,y,1);





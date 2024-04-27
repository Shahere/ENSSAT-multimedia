% MP3 = mp3read2

function TP_PERROT_BARBOTAUD
    function draw_spectrogramm(data, Fs)
        sigIn = data(1:end,1);
        figure
        subplot(211);
        t=(0:length(sigIn)-1)*(1/Fs);
        plot(t,sigIn);
        axis([t(1) t(end) min(sigIn) max(sigIn)]);
        xlabel('Discrete time')
        ylabel('Amplitude')
        subplot(212);
        L=floor(30e-3*Fs);
        Nfft=2^nextpow2(L);
        R = L-20;
        [sp,f] = spectrogram(sigIn, L, R, Nfft, Fs);
        t=(0:size(sp,2)-1) * ((L-R)/Fs);
        imagesc(t, f, 10*log10(abs(sp).^2))
        axis xy;
        xlabel('Discrete time')
        ylabel('Frequence (Hz)')
    end
    function spectrogram_song(data,Fs)
        sigIn = data(1:end,1);
        segment = sigIn(95697:111970);
        Nfft = 4*length(segment);
        spectre = fft(segment, Nfft).^2;
        f=(0:Nfft-1)*(Fs/Nfft);
        figure
        subplot(211);
        plot(segment)
        axis([1 length(segment) min(segment) max(segment)])
        xlabel('Discrete Time');ylabel('Amplitude')
        title('Représentation temporelle du segment audio')

        subplot(212);
        plot(f(1:1+Nfft/2),10*log10(abs(spectre(1:1+Nfft/2))));
        xlabel('Frequence (Hz)');ylabel('DSP (dB)');
        axis([0 Fs/2 -60 45])
        title('Représentation fréquentielle du segment audio')
    end
    function filtre(data,Fs)
        a11 = 0.001;
        a12 = -(1-a11);
        visuFiltre(a11,[1 -a12],Fs);
    end
    function y = applyFilter(data)
        sigIn = data(1:end,1);
        a11 = 0.001;
        a12 = 1-a11;
        x = sigIn;
        Zin = [];
        for k=1:length(x)
            [Zout, y(k)] = filtrage(a11, [1 -a12], x(k), Zin);
            Zin = Zout;
        end;
    end
    function wawa(data, Fs)
        sigIn = data(1:end,1);
        ym = data(1:end,1);
        yl = data(1:end,1);
        prec1 = 0;
        prec2 = 0;
        dump = 0.1;
        beta12 = -2*dump;
        beta13 = -1;
        Zin1 = eps;
        Zin2 = eps;
        fmin = 2000;
        fmax = 5000;
        fw = 3000;
        deltaf = fw / Fs;
        sigOut = data(1:end,1);

        alpha12 = 1;
        alpha22 = 1;
        fc = fmin: deltaf :fmax;
        while (length(fc) < length(sigIn))
            fc = [fc (fmax:-deltaf:fmin)];
            fc = [fc (fmin:+deltaf:fmax)];
        end
        alpha11 = 2 * sin(pi * fc / Fs);
        alpha21 = alpha11;

        ym(1) = alpha11(1) * sigIn(1);
        yl(1) = alpha21(1) * ym(1);

        for n=1:length(sigIn)
            yh(n)= sigIn(n)+beta12*prec1+beta13*prec2;

            [Zout1,ym(n)]= filtrage(alpha11(n),[1 -alpha12],yh(n),Zin1);
            [Zout2,yl(n)]= filtrage(alpha21(n),[1 -alpha22],ym(n),Zin2);
            Zin1 = Zout1;
            Zin2 = Zout2;
            prec1 = ym(n);
            prec2 = yl(n);
            sigOut(n) = yh(n) + ym(n) + yl(n);
        end
        soundsc(ym, Fs);
    end
%[data,Fs]=audioread('in_guitare.mp3');
%wawa(data, Fs)
%filtre(data, Fs)
%y = applyFilter(data);
%draw_spectrogramm(data,Fs);
%spectrogram_song(data,Fs)
%filtre(data, Fs)
% Le filtrage passe-bas laisse passer uniquement les basses fréqunces
end


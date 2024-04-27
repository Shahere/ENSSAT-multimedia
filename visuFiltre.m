function visuFiltre(Num,Den,fe)

[H,W] = freqz(Num,Den,2048,fe);
figure;
plot(W,20*log10(abs(H)))
grid;
axis([W(1) W(end) min(20*log10(abs(H))) max(20*log10(abs(H)))])
xlabel('Fréquence (Hz)')
ylabel('Atténuation (dB)');
title('Fonction de transfert du filtre')


function [Y,t] = complex2real(F,t)
%Y = complex2real(F,[t])
%
%Returns the real-valued amplitudes and phases in a structure calcullated
%from the complex-valued vector F having the convention of fft's output.
%This is the inverse of the function 'real2complex'.
%
%Inputs:
%   F        complex-valued vector in the convention of fft's output.
%   t        time vector of size y (default is 1:length(F));
%
%Outputs:    Structure Y with fields:
%   dc       mean value of y
%   amp      vector of amplitudes (length ceil(length(t)/2))
%   ph       vector of phases (in degrees, cosine phase)
%   nt       length of t (needed for myifft)
%
%SEE ALSO    real2complex fft ifft
%
%Example:
%
%t = 0:.01:.99;
%y=t<.5; 
%Y = complex2real(fft(y),t);
%clf;subplot(1,2,1);stem(t,y);
%xlabel('Time (s)');
%subplot(1,2,2);stem(Y.freq,Y.amp);
%xlabel('Frequency (Hz)')

%4/15/09     Written by G.M. Boynton at the University of Washington


%Deal with defaults
if ~exist('t','var')
    t = 1:length(F);
end

if isempty(t)
    t = 1:length(F);
end

%Calculate values based on t
nt = length(t);
dt = t(2)-t(1);

%DC is first value scaled by nt
dc = F(1)/nt;
%'real' amplitudes scale the fft by 2/nt
amp = 2*abs(F)/nt;
%'real' phases are reversed (and converted to degrees)
ph = -180*angle(F)/pi;

%Pull out the first half (omitting 'negative' frequencies)
id = 2:(ceil(nt/2)+1);
if size(F,2) == 1;
    id = id';
end

%Stuff the values in to the fields of Y
Y.dc = dc;
Y.ph = ph(id); %cosine phase
Y.amp = amp(id);
Y.freq = (1:length(id))/(nt*dt);
if size(F,2) == 1;
    Y.freq = Y.freq';
end

%Hack to deal with even vs. odd lengths of time series

if mod(nt,2) %length is odd
    Y.amp(end) = 0;
else %length is even
    Y.amp(end) = Y.amp(end)/2;
end

Y.nt= nt;




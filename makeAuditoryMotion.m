function [yy,p] = makeAuditoryMotion(p,s,x,y)
%[yy,p] = makeAuditoryMotion(p,s,x,y)
%
%Simulates sound signals for left and right ear given a monaural sound signal 
%and a vector of spatial position over time.  The sound intensity at each ear 
%is estimated by interpolating the sound signal with the time-varying delay
%according to the source's distance and the speed of sound.
%
%This automatically incorporates interaural time delay, inverse square law 
%attenuation, and doppler shift. Interaural level difference is partially 
%implemented by the inverse square law, but without a head shadow.
%
%Inputs:
% structure 'p' with fields:
%  required:
%   p.Fs    sampling rate for sound signal  (default is 44000 kHz)
%   s       monaural sound signal (column vector)
%  optional:
%   p.c     speed of sound (meters/second)  (default is 345)
%   p.a     width of head (meters)          (default is 0.0875)
%   p.dur   duration of sound (seconds)     (calculated if not provided)
%
% column vector 's' sound signal
% column vectors 'x' and 'y' of length(s)containing location of sound over time
%   in meters, with x=0, y=0 the location of center of head. 
%   x<0 left, and y>0 in front of head
%
%Output:
%  yy:  length(s) x 2 sound file for left and right channels
%  structure p with default parameters filled in.
%
%Written by G.M. Boynton, Spring of 2008 at the University of Washington

%calculate/generate default parameters if needed
if ~isfield(p,'c')
    p.c = 345;  %speed of sound (meters/second)
end

if ~isfield(p,'a')
    p.a = .0875;  %with of head (meters)
end

if ~isfield(p,'Fs')
    p.Fs = 44000;  %sampling rate for sound file (Hz)
end

if ~isfield(p,'dur')
    p.dur = p.Fs*length(s); %duration of sound file (seconds)
end

%time vector
t= linspace(0,p.dur,p.dur*p.Fs)';

%distance (squared) of sound source from left and right ears
if ne(size(x,2),size(y,2))
    y=y';
end;
size(x)
size(y)
dLsq = (x-p.a).^2+y.^2;
dRsq = (x+p.a).^2+y.^2;

%time delay (in seconds) from source to left and right ears
tauL = sqrt(dLsq)/p.c;  
tauR = sqrt(dRsq)/p.c;

%zero out the auditory stimulus matrix
yy = zeros(length(t),2);

%Here's the interesting part: 
%Interpolate back in time to determine sound for each ear and
%attenuate the sound based on inverse-square law
if ne(size(t,2),size(tauR,2))
    tauR=tauR';
end;
if ne(size(t,2),size(tauL,2))
    tauL=tauL';
end;
if ne(size(t,2),size(dRsq,2))
    dRsq=dRsq';
end;
if ne(size(t,2),size(dLsq,2))
    dLsq=dLsq';
end;
if ne(size(t,2),size(s,2))
    s=s';
end;
yy(:,1) = interp1(t,s./dRsq,t-tauR);
yy(:,2) = interp1(t,s./dLsq,t-tauL);

%clear out the NaN's (at the beginning of time where sound hasn't come to
%the head yet)
yy(isnan(yy))= 0;






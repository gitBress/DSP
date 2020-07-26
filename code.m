clear all 
close all
clc

[data,fs] = audioread('signal_028.wav');

fprintf('1 Original signal with Sampling Frequency = %d Hz\n',fs);
%%

%Time vector
signal_length = length(data)/fs;
t = linspace(0, signal_length, length(data));

%DFT of the signal
fprintf('2 Computation of the DFT of the original signal\n');

X = fft(data);
f = linspace(0, fs, length(X));

figure
subplot(311)
plot(t, data)
title('Raw signal')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([0 signal_length])
subplot(312)
plot(f,abs(X))
title('Magnitude of the DFT of the signal')
xlabel('Frequency (Hz)')
ylabel('|X(f)|')
xlim([0 fs])
subplot(313)
plot(f,mag2db(abs(X)))
title('Magnitude of the DFT of the signal in dB')
xlabel('Frequency (Hz)')
ylabel('|X(f)|dB')
xlim([0 fs])

%%
%Find peaks
MPD = 17000;
DIR = 'descend';
NP = 2;

[PKS,LOCS] = findpeaks(abs(X(1:(length(f)/2))),'MinPeakDistance',MPD,'SortStr',DIR,'NPeaks',NP);

%Plot peaks
figure
findpeaks(abs(X(1:(length(f)/2))),'MinPeakDistance',MPD,'SortStr',DIR,'NPeaks',NP);
title('Peaks found in the DFT of the signal')

%%
%Frequencies of sinusoidal carriers
f1 = (LOCS(1)*fs)/length(f);
f2 = (LOCS(2)*fs)/length(f);

fprintf('3 First peak found at frequency = %5.0f Hz\n',f1);
fprintf('  Second peak found at frequency = %5.0f Hz\n',f2);

%Estimate of the amplitudes of the carriers
A1=2*abs(PKS(1))/length(f);
A2=2*abs(PKS(2))/length(f);

fprintf('4 Estimated amplitude of the first carrier = %1.2f\n',A1);
fprintf('  Estimated amplitude of the second carrier = %1.2f\n',A2);
%%
fprintf('5 Extraction of the two carriers with second order IIR filters\n');

delta_f_3db = 10; %Hz

%Extraction of the first carrier at f1=14.3kHz
[Num1 Den1] = second_order_IIR(delta_f_3db, fs, f1);
carrier1 = filter(Num1, Den1, data);
figure
freqz(Num1, Den1, length(f), 'whole', fs);
title('Frequency response of the filter to extract the carrier at 14.3kHz')

fprintf('\tSPECIFICATIONS OF THE 1st SECOND ORDER IIR FILTER\n');
fprintf('\t3dB Bandwidth = %2.0f Hz\n',delta_f_3db);
fprintf('\tCenter frequency = %5.0f Hz\n',f1);


%Extraction of the second carrier at f2=31.8kHz
[Num2 Den2] = second_order_IIR(delta_f_3db, fs, f2);
carrier2 = filter(Num2, Den2, data);
figure
freqz(Num2, Den2, length(f), 'whole', fs);
title('Frequency response of the filter to extract the carrier at 31.8kHz')

fprintf('\tSPECIFICATIONS OF THE 2nd SECOND ORDER IIR FILTER\n');
fprintf('\t3dB Bandwidth = %2.0f Hz\n',delta_f_3db);
fprintf('\tCenter frequency = %5.0f Hz\n',f2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(211)
% plot(t,carrier1)
% title('First carrier extracted, f1=14.3kHz')
% xlim([0 signal_length])
% ylim([-0.1 0.1])
% xlabel('Time [s]')
% ylabel('Amplitude')
% subplot(212)
% plot(t,carrier2)
% title('Second carrier extracted, f2=31.8kHz')
% xlim([0 signal_length])
% ylim([-0.1 0.1])
% xlabel('Time [s]')
% ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%demodulation of the two information signals
x1 = (1/A1).*(carrier1.*data);
x2 = (1/A2).*(carrier2.*data);
%(multiplication by a constant to get the output audible)

fprintf('6 Demodulation of the two information signals\n');
%%
fprintf('7 Filtering of the demodulated signal to reject audible distortions\n');

%Filtering demodulated signals: LOW-PASS

%Specifications of the low-pass filter
rp = 3;           % Passband ripple [dB]
rs = 80;          % Attenuation at the stopband [dB]
fp = [8000 9000]; % Cutoff frequencies
a = [1 0];       % Desired amplitudes

fprintf('\tSPECIFICATIONS OF THE LOW-PASS LINEAR-PHASE FIR FILTER\n');
fprintf('\tPass-band frequency = %4.0f Hz\n',fp(1));
fprintf('\tStop-band frequency = %4.0f Hz\n',fp(2));
fprintf('\tDesired amplitude of pass-band = %1.0f \n',a(1));
fprintf('\tDesired amplitude of stop-band = %1.0f \n',a(2));
fprintf('\tAttenuation at the stop-band = %1.0f dB\n',rs);
fprintf('\tPass-band ripple = %1.0f dB\n',rp);

%Parameters computation
fo = [0 2*fp(1)/fs 2*fp(2)/fs 1];
ao = [a(1) a(1) a(2) a(2)];
n = -(10*log(10^(rp/20)*10^(-rs/20))+13)/(14.6*(fp(2)-fp(1))/fs); %order of the filter
%Filter
b = firpm(ceil(n),fo,ao);
figure
freqz(b,1,length(f),'whole',fs);
title('Frequency response of the low-pass filter')


%Low-pass filtering
x1_lowpass =  filter(b, 1, x1);
x2_lowpass =  filter(b, 1, x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(221)
% plot(t,x1)
% title('First demodulated signal')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
% subplot(222)
% plot(t,x1_lowpass)
% title('First demodulated signal after low-pass filtering')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
% 
% subplot(223)
% plot(t,x2)
% title('Second demodulated signal')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
% subplot(224)
% plot(t,x2_lowpass)
% title('Second demodulated signal after low-pass filtering')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Filtering demodulated signals: HIGH-PASS NOTCH IIR

%Specifications of the high-pass filter
f0 = 0; %Hz
delta_f_3db_notch = 40; %Hz

fprintf('\tSPECIFICATIONS OF THE HIGH-PASS NOTCH IIR FILTER\n');
fprintf('\tCenter frequency = %1.0f Hz\n',f0);
fprintf('\t3dB Bandwidth = %2.0f Hz\n',delta_f_3db_notch);

% Computation of the notch filter coefficients
b1 = -2*cos(2*pi*f0/fs); 
Num_notch = [1 b1 1];
delta = pi*delta_f_3db_notch/fs; r = 1-delta;
a2 = r^2; a1 =-2*r*cos(2*pi*f0/fs);
Den_notch= [1 a1 a2];

freqz(Num_notch,Den_notch,length(f),'whole',fs);
title('Frequency response of the high-pass filter')

%High-pass filtering
x1_demod =  filter(Num_notch,Den_notch, x1_lowpass);
x2_demod =  filter(Num_notch,Den_notch, x2_lowpass);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(221)
% plot(t,x1_lowpass)
% title('First demodulated signal low-pass filtered')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
% subplot(222)
% plot(t,x1_demod)
% title('First demodulated signal band-pass filtered')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
% 
% subplot(223)
% plot(t,x2_lowpass)
% title('Second demodulated signal low-pass filtered')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
% subplot(224)
% plot(t,x2_demod)
% title('Second demodulated signal band-pass filtered')
% xlim([0 signal_length])
% xlabel('Time [s]')
% ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Spectrum of the two demodulated signals
fprintf('8 Computation of the spectrum of the two demodulated signals\n');

Y1 = fft(x1_demod);
Y2 = fft(x2_demod);

figure
subplot(221)
plot(f,abs(Y1))
title('DFT of the first demodulated signal')
xlabel('Frequency (Hz)')
ylabel('|Y1(f)|')
xlim([0 fs])
subplot(223)
plot(f,mag2db(abs(Y1)))
title('DFT of the first demodulated signal in dB')
xlabel('Frequency (Hz)')
ylabel('|Y1(f)|dB')
xlim([0 fs]) 
subplot(222)
plot(f,abs(Y2))
title('DFT of the second demodulated signal')
xlabel('Frequency (Hz)')
ylabel('|Y2(f)|')
xlim([0 fs]) 
subplot(224)
plot(f,mag2db(abs(Y2)))
title('DFT of the second demodulated signal in dB')
xlabel('Frequency (Hz)')
ylabel('|Y2(f)|dB')
xlim([0 fs]) 


%%
fprintf('9 Writing of the stereo audio file containing the two output signals\n');

audiowrite('bressan_giulia.wav', [x1_demod x2_demod], fs);

%%
%Filter design ---> Second Order IIR filter
function [B A] = second_order_IIR(band, fs, f0)
%fixed order of the filter
delta = (pi*band)/(fs);
r = 1-delta;
%b_0 = 2*(1-r)*abs(sin(2*pi*f0/fs));
b_0 = (1-r); %with two zeros at numerator
%coefficients of the transfer function
A = [1 -2*cos(2*pi*f0/fs) r^2]; %denominator of the transfer function
B = [b_0 0 -b_0]; %numerator of the transfer function
end

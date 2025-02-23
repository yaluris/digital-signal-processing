clear;
close all;

%% 1
Wc = 0.4*pi;
Fc = Wc/(2*pi);
Fs = 100;
N = 21;

Wn = Fc/(Fs/2);                             %normalized cutoff frequency

rec_filter = fir1(N-1, Wn, rectwin(N));     %rectangular filter

ham_filter = fir1(N-1, Wn, hamming(N));     %hamming filter

%frequency response
[h1, w1] = freqz(rec_filter, N);
[h2, w2] = freqz(ham_filter, N);

figure;
plot(w1, abs(h1), w2, abs(h2));
legend('Rectangular', 'Hamming');
xlabel('Frequency');
ylabel('Magnitude');
title('Frequency Response');

%% 2A
Wc = 0.5*pi;
Fs = 100;
Fc = Wc/(2*pi);

Wn = Fc/(Fs/2);                                %normalized cutoff freq

Na = 21; 
Nb = 41;

ham1_filter = fir1(Na-1, Wn, hamming(Na));     %hamming filter
ham2_filter = fir1(Nb-1, Wn, hamming(Nb));     %hamming filter

[ha, wa] = freqz(ham1_filter, Na);
[hb, wb] = freqz(ham2_filter, Nb);

figure;
subplot(1,2,1);
plot(wa, abs(ha));
xlabel('Frequency');
ylabel('Magnitude');
title('Hamming Frequency Response (N=21)');

subplot(1,2,2);
plot(wb, abs(hb));
xlabel('Frequency');
ylabel('Magnitude');
title('Hamming Frequency Response (N=41)');

han1_filter = fir1(Na-1, Wn, hann(Na));         %hanning filter
han2_filter = fir1(Nb-1, Wn, hann(Nb));         %hanning filter

[hc, wc] = freqz(han1_filter, Na);
[hd, wd] = freqz(han2_filter, Nb);

figure;
subplot(1,2,1);
plot(wc, abs(hc));
xlabel('Frequency');
ylabel('Magnitude');
title('Hanning Frequency Response (N=21)');

subplot(1,2,2);
plot(wd, abs(hd));
xlabel('Frequency');
ylabel('Magnitude');
title('Hanning Frequency Response (N=41)');

%% 2B
%Create signal x(t) = sin(15t)+0.25sin(200t)
figure;
t = [0:0.0001:5];
xt = sin(15*t)+0.25*sin(200*t);
plot(t,xt);
hold on;
Fs = 100; %Sampling frequency
Ts = 1/Fs; %Sampling period
N = 500; %Total samples are 500
n = 0:N-1;
x = sin(15*n*Ts)+0.25*sin(200*n*Ts); %Sampled x(t)
stem(n*Ts,x);
xlabel('t');
ylabel('x(t)');
title('Sampling x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz');
grid on;
hold off;
%According to the Nyquist sampling theorem, a signal can be perfectly reconstructed from its samples if the sampling frequency is greater than or equal to twice the maximum frequency
%x(t) = sin(15t)+0.25sin(200t) = sin(2*pi*2.39t)+0.25sin(2*pi*31.83t)
%Here Fs = 100Hz >= 2Fmax = 2*31.83 = 63.66Hz (2Fmax: Nyquist frequency)
%So there is no aliasing

%Signal spectrum (fasma) is the graphical representation of the amplitude of the Fourier transform for 1 period, which is around 0
%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz
figure;
F = [-Fs/2:Fs/N:Fs/2-Fs/N]; %Frequency vector (1 period => Fs/2, minus 1 sample because of 0 => Fs/2-Fs/N, divide the samples into Fs frequencies => step Fs/N)
Xf = fftshift(fft(x)); %Fourier transform x and rearrange it by shifting the zero-frequency component to the center (so that the entire spectrum appears)
plot(F,abs(Xf));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz');
grid on;

%Hamming filter with length N=21
X_hamfilter1 = filter(ham1_filter,1,x);
Xf_hamfilter1 = fftshift(fft(X_hamfilter1));
%Hamming filter with length N=41
X_hamfilter2 = filter(ham2_filter,1,x);
Xf_hamfilter2 = fftshift(fft(X_hamfilter2));
%Hanning filter with length N=21
X_hanfilter1 = filter(han1_filter,1,x);
Xf_hanfilter1 = fftshift(fft(X_hanfilter1));
%Hanning filter with length N=41
X_hanfilter2 = filter(han2_filter,1,x);
Xf_hanfilter2 = fftshift(fft(X_hanfilter2));

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hamming filter with length N=21
figure;
subplot(2,1,1);
plot(F,abs(Xf_hamfilter1));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hamfilter1|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz after applying Hamming filter with length N=21');
grid on;

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hamming filter with length N=41
subplot(2,1,2);
plot(F,abs(Xf_hamfilter2));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hamfilter2|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz after applying Hamming filter with length N=41');
grid on;

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hanning filter with length N=21
figure;
subplot(2,1,1);
plot(F,abs(Xf_hanfilter1));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hanfilter1|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz after applying Hanning filter with length N=21');
grid on;

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hanning filter with length N=41
subplot(2,1,2);
plot(F,abs(Xf_hanfilter2));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hanfilter2|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 100Hz after applying Hanning filter with length N=41');
grid on;

%% 2C
figure;
t = [0:0.0001:10];
xt = sin(15*t)+0.25*sin(200*t);
plot(t,xt);
hold on;
Fs = 50;
Ts = 1/Fs;
N = 500;
n = 0:N-1;
x = sin(15*n*Ts)+0.25*sin(200*n*Ts);
stem(n*Ts,x);
xlabel('t');
ylabel('x(t)');
title('Sampling x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz');
grid on;
hold off;
%x(t) = sin(15t)+0.25sin(200t) = sin(2*pi*2.39t)+0.25sin(2*pi*31.83t)
%Here Fs = 50Hz < 2Fmax = 2*31.83 = 63.66Hz (2Fmax: Nyquist frequency)
%So there is aliasing

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz
figure;
F = [-Fs/2:Fs/N:Fs/2-Fs/N];
Xf = fftshift(fft(x));
plot(F,abs(Xf));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz');
grid on;

%Hamming filter with length N=21
X_hamfilter1 = filter(ham1_filter,1,x);
Xf_hamfilter1 = fftshift(fft(X_hamfilter1));
%Hamming filter with length N=41
X_hamfilter2 = filter(ham2_filter,1,x);
Xf_hamfilter2 = fftshift(fft(X_hamfilter2));
%Hanning filter with length N=21
X_hanfilter1 = filter(han1_filter,1,x);
Xf_hanfilter1 = fftshift(fft(X_hanfilter1));
%Hanning filter with length N=41
X_hanfilter2 = filter(han2_filter,1,x);
Xf_hanfilter2 = fftshift(fft(X_hanfilter2));

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hamming filter with length N=21
figure;
subplot(2,1,1);
plot(F,abs(Xf_hamfilter1));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hamfilter1|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz after applying Hamming filter with length N=21');
grid on;

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hamming filter with length N=41
subplot(2,1,2);
plot(F,abs(Xf_hamfilter2));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hamfilter2|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz after applying Hamming filter with length N=41');
grid on;

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hanning filter with length N=21
figure;
subplot(2,1,1);
plot(F,abs(Xf_hanfilter1));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hanfilter1|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz after applying Hanning filter with length N=21');
grid on;

%Spectrum of signal x(t) = sin(15t)+0.25sin(200t) after applying Hanning filter with length N=41
subplot(2,1,2);
plot(F,abs(Xf_hanfilter2));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hanfilter2|');
title('Spectrum of signal x(t) = sin(15t)+0.25sin(200t) with Fs = 50Hz after applying Hanning filter with length N=41');
grid on;

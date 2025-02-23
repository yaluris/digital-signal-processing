clear;
close all;

%% 1
%lowpass Butterworth filter
fs = 10000;                         %sampling freq
            
Wp = 2*pi*3000;                     %normalized passband freq
Ws = 2*pi*4000;                     %normalized stopband freq   

[n,Wn] = buttord(Wp,Ws,3,30,'s');   %n order of analog filter

%poles and gain of n order Butterworth filter
[z,p,k] = buttap(n);    

%convert zero-pole-gain filter parameters to transfer function form
[b,a] = zp2tf(z,p,k);

%change cutoff frequency for lowpass analog filter
[bt,at] = lp2lp(b,a,Wn);

N = 2048;                           %number of samples    
f = 0:fs/(N-1):fs/2;                %frequency axis
w = 2*pi*f;                         %angular frequencies

%frequency response of analog filters
analog_f = freqs(bt,at,w);

%for analog-to-digital filter conversion
[num1, den1] = bilinear(bt, at, fs);

%frequency response of digital filter
digital_f = freqz(num1, den1, f, fs);

figure;
plot(f, mag2db(abs(analog_f)), '--', f, mag2db(abs(digital_f)));
legend('Analog Filter', 'Digital Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Butterworth Lowpass Filter');
ylim([-350 30]);

%attenuation stopband 50db 
%Rs = 50; 
[n,Wn] = buttord(Wp,Ws,3,50,'s');   %n order of analog filter

%poles and gain of n order Butterworth filter
[z,p,k] = buttap(n);    

%convert zero-pole-gain filter parameters to transfer function form
[b,a] = zp2tf(z,p,k);

%change cutoff frequency for lowpass analog filter
[bt,at] = lp2lp(b,a,Wn);

N = 2048;                           %number of samples    
f = 0:fs/(N-1):fs/2;                %frequency axis
w = 2*pi*f;                         %angular frequencies

%frequency response of analog filters
analog_f = freqs(bt,at,w);

%for analog-to-digital filter conversion
[num2, den2] = bilinear(bt, at, fs);

%frequency response of digital filter
digital_f = freqz(num2, den2, f, fs);

figure;
plot(f, mag2db(abs(analog_f)), '--', f, mag2db(abs(digital_f)));
legend('Analog Filter', 'Digital Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Butterworth Lowpass Filter');
ylim([-350 30]);

%% 2
Wc = 2;                         %cutoff (rad/sec)
Fc = Wc/(2*pi);                 %cutoff freq (Hz)
Ts = 0.2;                       %sampling period (sec)
Fs = 1/Ts;                      %sampling freq (Hz)
F_nyquist = Fs/2;               %nyquist freq
pb_r = 3 ;                      %passband ripple (dB)
F_cheby1 = Fc/F_nyquist;        %passband edge frequency

N = 256;                        %samples
Fa = 0: 1/(N-1) :1;             %normalized frequency axis

ftype = 'high';                 %specifies a highpass filter with passband edge frequency
n1 = 2;                         %order 2
n2 = 16;                        %order 16

[b, a] = cheby1(n1, pb_r, F_cheby1 ,ftype);     %designs a Chebyshev Type I filter, depending on the value of ftype and the number of elements of F_cheby1(here highpass, and Fc/F_nyquist respectively ).
highpass_filter_of_n1_order = freqz(b, a, N);   %returns the n-point frequency response vector h and the corresponding angular frequency vector w for the digital filter with transfer function coefficients stored in b and a.

[d, c] = cheby1(n2, pb_r, F_cheby1 ,ftype);     %designs a Chebyshev Type I filter, depending on the value of ftype and the number of elements of F_cheby1(here highpass, and Fc/F_nyquist respectively ).
highpass_filter_of_n2_order = freqz(d, c, N);   %returns the n-point frequency response vector h and the corresponding angular frequency vector w for the digital filter with transfer function coefficients stored in b and a.

figure;
plot(Fa, mag2db(abs(highpass_filter_of_n1_order)), ':', Fa, mag2db(abs(highpass_filter_of_n2_order))); %expresses in decibels (dB) the magnitude measurements specified in highpass_filter_of_n1_order. The relationship between magnitude and decibels is ydb = 20 log10(highpass_filter_of_n1_order).
title('Highpass Chebyshev Type I Filter');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
legend('n=2', 'n=16');

%% 3A
%Create signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t)
figure;
t = [0:0.0001:0.05];
xt = 1+cos(1000*t)+cos(16000*t)+cos(30000*t);
plot(t,xt);
hold on;
Fs = 10000; %Sampling frequency
Ts = 1/Fs; %Sampling period
N = 500; %Total samples are 500
n = 0:N-1;
x = 1+cos(1000*n*Ts)+cos(16000*n*Ts)+cos(30000*n*Ts); %Sampled x(t)
stem(n*Ts,x);
xlabel('t');
ylabel('x(t)');
title('Sampling x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t)');
grid on;
hold off;
%According to the Nyquist sampling theorem, a signal can be perfectly reconstructed from its samples if the sampling frequency is greater than or equal to twice the maximum frequency
%x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t) = 1+cos(2*pi*159t)+cos(2*pi*2547t)+cos(2*pi*4775t)
%Here Fs = 10000Hz >= 2Fmax = 2*4775 = 9550Hz (2Fmax: Nyquist frequency)
%So there is no aliasing

%Signal spectrum (fasma) is the graphical representation of the amplitude of the Fourier transform for 1 period, which is around 0
%Spectrum of signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t)
figure;
F = [-Fs/2:Fs/N:Fs/2-Fs/N]; %Frequency vector (1 period => Fs/2, minus 1 sample because of 0 => Fs/2-Fs/N, divide the samples into Fs frequencies => step Fs/N)
Xf = fftshift(fft(x)); %Fourier transform x(t) and rearrange it by shifting the zero-frequency component to the center (so that the entire spectrum appears)
plot(F,abs(Xf));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf|');
title('Spectrum of signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t)');
grid on;

%Spectrum of signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t) after applying lowpass Butterworth filter with stopband attenuation 30db
figure;
X_lpfilter1 = filter(num1,den1,x);
Xf_lpfilter1 = fftshift(fft(X_lpfilter1));
plot(F,abs(Xf_lpfilter1));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_lpfilter1|');
title('Spectrum of signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t) after applying lowpass Butterworth filter with stopband attenuation 30db');
grid on;

%Spectrum of signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t) after applying lowpass Butterworth filter with stopband attenuation 50db
figure;
X_lpfilter2 = filter(num2,den2,x);
Xf_lpfilter2 = fftshift(fft(X_lpfilter2));
plot(F,abs(Xf_lpfilter2));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_lpfilter2|');
title('Spectrum of signal x(t) = 1+cos(1000t)+cos(16000t)+cos(30000t) after applying lowpass Butterworth filter with stopband attenuation 50db');
grid on;

%% 3B
%Create signal x(t) = 1+cos(1.5t)+cos(5t)
figure;
t = [0:0.0001:100];
xt = 1+cos(1.5*t)+cos(5*t);
plot(t,xt);
hold on;
Fs = 5; %Sampling frequency
Ts = 1/Fs; %Sampling period
N = 500; %Total samples are 500
n = 0:N-1;
x = 1+cos(1.5*n*Ts)+cos(5*n*Ts); %Sampled x(t)
stem(n*Ts,x);
xlabel('t');
ylabel('x(t)');
title('Sampling x(t) = 1+cos(1.5t)+cos(5t)');
grid on;
hold off;
%x(t) = 1+cos(1.5t)+cos(5t) = 1+cos(2*pi*0.24t)+cos(2*pi*0.8t)
%Here Fs = 5Hz >= 2Fmax = 2*0.8 = 1.6Hz (2Fmax: Nyquist frequency)
%So there is no aliasing

%Spectrum of signal x(t) = 1+cos(1.5t)+cos(5t)
figure;
F = [-Fs/2:Fs/N:Fs/2-Fs/N]; %Frequency vector (1 period => Fs/2, minus 1 sample because of 0 => Fs/2-Fs/N, divide the samples into Fs frequencies => step Fs/N)
Xf = fftshift(fft(x)); %Fourier transform x(t) and rearrange it by shifting the zero-frequency component to the center (so that the entire spectrum appears)
plot(F,abs(Xf));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf|');
title('Spectrum of signal x(t) = 1+cos(1.5t)+cos(5t)');
grid on;

%Spectrum of signal x(t) = 1+cos(1.5t)+cos(5t) after applying highpass Chebyshev filter with order 2
figure;
X_hpfilter1 = filter(b,a,x);
Xf_hpfilter1 = fftshift(fft(X_hpfilter1));
plot(F,abs(Xf_hpfilter1));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hpfilter1|');
title('Spectrum of signal x(t) = 1+cos(1.5t)+cos(5t) after applying highpass Chebyshev filter with order 2');
grid on;

%Spectrum of signal x(t) = 1+cos(1.5t)+cos(5t) after applying highpass Chebyshev filter with order 16
figure;
X_hpfilter2 = filter(d,c,x);
Xf_hpfilter2 = fftshift(fft(X_hpfilter2));
plot(F,abs(Xf_hpfilter2));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf_hpfilter2|');
title('Spectrum of signal x(t) = 1+cos(1.5t)+cos(5t) after applying highpass Chebyshev filter with order 16');
grid on;

clear;
close all;

%% 1A
%Create the signals x[n] = 4u[n-4]+r[-n+8] and y[n] = (1/2)^|n|
figure;
nx = [-15:15]; %Time of x[n]
X = 4*(nx-4>=0)+(-nx+8>=0).*(-nx+8);
subplot(2,1,1);
stem(nx,X);
xlabel('n');
ylabel('x[n]');
title('x[n] = 4u[n-4]+r[-n+8]');
grid on;

ny = [-15:10]; %Time of y[n]
Y = (1/2).^abs(ny);
subplot(2,1,2);
stem(ny,Y);
xlabel('n');
ylabel('y[n]');
title('y[n] = (1/2)^{|n|}');
grid on;

nh = [nx(1)+ny(1):nx(end)+ny(end)]; %Time of h[n] = x[n]*y[n] (convolution)
LenX = length(nx); %Length of x[n] is 31
LenY = length(ny); %Length of y[n] is 26
LenH = length(nh); %Length of convolution is LenX+LenY-1 = 56

%Zero padding X so that the signals are the same size and can be multiplied (fill X left and right with as many zeros as the length of Y minus 1 point, which is the minimum overlap)
figure;
X0 = [zeros(1,LenY-1) X zeros(1,LenY-1)];
k = [nx(1)-(LenY-1):nx(end)+(LenY-1)]; %Time of X0
subplot(2,1,1);
stem(k,X0);
xlabel('k');
ylabel('x[k]');
grid on;

Y_rev = Y(end:-1:1); %Anaklash y[-k]

%Discrete time convolution formula
for i = 1:LenH %Length of h[n]
    Y_rev0 = [zeros(1,(i-1)) Y_rev zeros(1,(LenH-i))]; %Zero padding Y_rev so that the signals are the same size and can be multiplied (the padding changes as the signal shifts to the right)
    subplot(2,1,2);
    stem(k,Y_rev0);
    xlabel('k');
    ylabel('y[nh(i)-k]');
    grid on;
    H1(i) = sum(X0.*Y_rev0); %The convolution h[n] = x[n]*y[n] for n = nh(i) is equal to the sum of the products of the overlapping elements
end

figure;
subplot(2,1,1);
stem(nh,H1);
xlabel('n');
ylabel('h[n]');
title('Convolution without using the conv function');
grid on;

H2 = conv(X,Y);
subplot(2,1,2);
stem(nh,H2);
xlabel('n');
ylabel('h[n]');
title('Convolution using the conv function');
grid on;

%% 1B
n1 = -6:6;

xn = 1/4.*exp(abs(n1));
hn = 1/4*(n1>=-6);

figure;
stem(n1,xn)
title('xn = 1/4.*exp(abs(n1))');
xlabel('n1');
ylabel('x[n]');
figure;
stem(n1,hn)
title('hn = 1/4*(n1>=-6)');
xlabel('n1');
ylabel('h[n]');

yn = conv(xn,hn);
n = 2*n1(1):2*n1(end);
figure;
stem(n,yn)
title('Convolution');
xlabel('n');
ylabel('y[n]');


f = length(n);
xf = fft(xn, f);
hf = fft(hn, f);

multf = xf.*hf;

imf = ifft(multf); %ifft, return to time-domain diagram
figure;
stem(n, imf)
title('Convolution');
xlabel('n');
ylabel('imf[n]');

%% 2
t = 0 : 0.0001 : 0.5;
xt = 5*cos(24*pi*t) - 2*sin(1.5*pi*t);

figure;
plot(t, xt);
title('5*cos(24*pi*t) - 2*sin(1.5*pi*t))');
xlabel('t');
ylabel('x(t)');
grid on;

%Ts=1/48s
fs=48;
tSamp=0:1/fs:0.5;

xSamp = 5*cos(24*pi*tSamp) - 2*sin(1.5*pi*tSamp);
figure;
plot(t, xt);
title('Sampling Ts=1/48s');
xlabel('Time(s)');
hold on
grid on
plot(tSamp,xSamp)
%stem(tSamp,xSamp)
legend('x(t)', 'sampling of x(t)')
hold off

%Ts=1/24s
fs2=24;
tSamp2=0:1/fs2:0.5;
xSamp2= 5*cos(24*pi*tSamp2) - 2*sin(1.5*pi*tSamp2);
figure;
plot(t, xt);
title('Sampling Ts=1/24s');
xlabel('Time(s)');
hold on
grid on
plot(tSamp2,xSamp2)
%stem(tSamp2,xSamp2)
legend('x(t)', 'sampling of x(t)')
hold off

%Ts=12s
fs3=12;
tSamp3=0:1/fs3:0.5;
xSamp3= 5*cos(24*pi*tSamp3) - 2*sin(1.5*pi*tSamp3);
figure;
plot(t, xt);
title('Sampling Ts=1/12s');
xlabel('Time(s)');
hold on
grid on
plot(tSamp3,xSamp3)
%stem(tSamp3,xSamp3)
legend('x(t)', 'sampling of x(t)')
hold off

%Ts=1/46s
fs4=46;
tSamp4=0:1/fs4:0.5;
xSamp4= 5*cos(24*pi*tSamp4) - 2*sin(1.5*pi*tSamp4);
figure;
plot(t, xt);
title('Sampling Ts=1/46s');
xlabel('Time(s)');
hold on
grid on
plot(tSamp4,xSamp4)
%stem(tSamp4,xSamp4)
legend('x(t)', 'sampling of x(t)')
hold off

%% 3A
%Create signal x(t) = 10cos(2*pi*20t)-4sin(2*pi*40t+5)
figure;
t = [0:0.0001:0.64];
xt = 10*cos(2*pi*20*t)-4*sin(2*pi*40*t+5);
plot(t,xt);
hold on;
%According to the Nyquist sampling theorem, a signal can be perfectly reconstructed from its samples if the sampling frequency is greater than or equal to twice the maximum frequency
%So to avoid aliasing Fs >= 2Fmax => Fs >= 2*40 = 80Hz (2Fmax: Nyquist frequency)
Fs = 200; %Sampling frequency
Ts = 1/Fs; %Sampling period
N = 128; %Total samples are 128
n = 0:N-1;
x = 10*cos(2*pi*20*n*Ts)-4*sin(2*pi*40*n*Ts+5); %Sampled x(t)
stem(n*Ts,x);
xlabel('t');
ylabel('x(t)');
title('Sampling x(t) = 10cos(2*pi*20t)-4sin(2*pi*40t+5)');
grid on;
hold off;

%Signal spectrum (fasma) is the graphical representation of the amplitude of the Fourier transform for 1 period, which is around 0
%Spectrum of signal x(t) = 10cos(2*pi*20t)-4sin(2*pi*40t+5)
figure;
F = [-Fs/2:Fs/N:Fs/2-Fs/N]; %Frequency vector (1 period => Fs/2, minus 1 sample because of 0 => Fs/2-Fs/N, divide the samples into Fs frequencies => step Fs/N)
Xf = fftshift(fft(x)); %Fourier transform x(t) and rearrange it by shifting the zero-frequency component to the center (so that the entire spectrum appears)
plot(F,abs(Xf));
xlabel('F(Hz)');
ylabel('Fourier Amplitude |Xf|');
title('Spectrum of signal x(t) = 10cos(2*pi*20t)-4sin(2*pi*40t+5)');
grid on;

%% 3B
%100:125:475 (from:step:to)
fs = 8000;
Ts = 1/fs;
theta = 300;
t = -0.01: Ts : 0.01-Ts;                        %Vector of time
figure;

for f0 = 100:125:475
    xt = sin(2*pi*f0*t + theta);
    NFFT = 2^nextpow2(length(xt));               %Number of Fourier samples
    XF = fftshift(fft(xt, NFFT)*Ts);             %Normalization of Fourier, because  we want to approach the integral, with sum
    freq = (-fs/2) : fs/NFFT : (fs/2- fs/NFFT) ; %Plot aplitude spectrum
    subplot(2,2,(fix(f0/100)));
    plot( freq, abs(XF) )                        %We use abs for keeping only the positive amplitudes of fourier
    if((f0/100) == 1)
        title(  '                             X(F) = (i*[(Pi/2)^(1/2)]*e^(-i*theta))*[ δ[f -2*pi*f0*Ts] ) - e^(2*i*theta)*δ[f + 2*pi*f0*Ts]) ]' );
    end
    xlabel( 'Frequency(Hz)' );
    ylabel( 'X(F)' );
end

% 7525:125:7900
figure;
k=1;

for f0 = 7525:125:7900
    xt = sin(2*pi*f0*t + theta);
    NFFT = 2^nextpow2(length(xt));               %Number of Fourier samples
    XF = fftshift(fft(xt, NFFT)*Ts);             %Normalization of Fourier, because  we want to approach the integral, with sum
    freq = (-fs/2) : fs/NFFT : (fs/2- fs/NFFT) ; %Plot aplitude spectrum
    subplot(2,2,k);
    plot( freq, abs(XF) )                        %We use abs for keeping only the positive amplitudes of fourier
    if(k == 1)
        title(  '                               X(F) = (i*[(Pi/2)^(1/2)]*e^(-i*theta))*[ δ[f - 2*pi*f0*Ts] ) - e^(2*i*theta)*δ[f + 2*pi*f0*Ts]) ]' );
    end  
    xlabel( 'Frequency(Hz)' );
    ylabel( 'X(F)' );
    k=k+1;
end

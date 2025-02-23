clear;
close all;

%% 1B
%H(z)=0.2z/(z^2-0.7z-0.18)
numerator=[0.2,0];
denominator=[1,-0.7,-0.18];
Hz=tf(numerator,denominator)
figure
zplane(numerator,denominator)
grid on
title('zero-pole');

%% 1D
n=-pi:pi/128:pi;
figure
freqz(numerator,denominator,n);
grid on
title('freqz(b,a,n)');

figure
freqz(numerator,denominator);
grid on
title('freqz(b,a)');

%% 1E
numerator2=[0.2,0];
denominator2=[1,-1.7,0.52,0.18];
Hz2=tf(numerator2,denominator2)
figure
freqz(numerator2,denominator2,n);
grid on
title('freqz(b,a,n)');

%% 2A
numerator = [4 -3.5 0];
denominator = [1 -2.5 1];

Hz = tf(numerator, denominator);                   %Hz = (4-(3.5.*(z^(-1))))/(1-(2.5.*(z^(-1)))+(z^(-2))), creates a continuous-time transfer function model, setting the numerator and denominator properties

[res, pol, ko] = residuez(numerator, denominator); %Finds the residues (A, B,..., N, direct terms of a partial fraction expansion of the ratio of numerator and denominator polynomials)
                                                   %Finds poles
                                                   %Direct terms, returned as a row vector. The direct term coefficient vector k is empty if length(b) is less than length(a). Otherwise: length(k) = length(b)-length(a)+1

A = res(1);                                        %We take the first number of our vector, which is our A for the partial fraction expansion of the ratio of numerator and denominator polynomials   
B = res(2);                                        %We take the second (and for this exercise last) number of our vector, which is our B for the partial fraction expansion of the ratio of numerator and denominator polynomials   
P1 = pol(1);                                       %We take the first number of our vector, which is our first pole
P2 = pol(2);                                       %We take the second (and for this exercise last) pole

syms z                                             %Lists the names of all symbolic scalar variables, functions, matrix variables, matrix functions and arrays in the workspace. Here, we declare a complex variable z to use it in the Z-Transform
T_Hz = A/(1- P1*z^(-1)) + B/(1- P2*z^(-1));
pretty(T_Hz)                                       %Prints T_Hz

%% 2B
inv_T_Hz = iztrans(T_Hz)

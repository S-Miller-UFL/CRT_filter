clear all; close all; clc;

n = 1:20;
w= -pi:pi/500:pi-pi/500;
impulse = zeros(1,20);
impulse(1) = 1;

%% filter using only poles

%% low pass filter(blur effect)
zeroeslow = [1];
poleslow = [1  -.9];
impulselow = filter(zeroeslow,poleslow,impulse);
%% high pass filter(distortion effect)
zeroeshigh = [1];
poleshigh = [1 .9];
impulsehigh = filter(zeroeshigh,poleshigh,impulse);
%% this high pass filter produces a really good distortion effect
%zeroes = [1];
%poles = [1 .9];
%impulse = filter(zeroes,poles,impulse);

%% filter using only zeroes
%% high pass filter(produces gray garbage)
%zeroes = [1 -1];
%poles = [1];
%impulse = filter(zeroes,poles,impulse);
%% low pass filter (does nothing noticeable)
%zeroes = [1 1];
%poles = [1];
%impulse = filter(zeroes,poles,impulse);
%figure("Name","magnitude response of input");
%plot(w,fft(x,length(w)));
%title("magnitude response of input")
%xlabel("frequency(radians)");
%ylabel("magnitude");
figure("Name","low pass plot")
pzplot(zeroeslow,poleslow);
figure("Name","high pass plot")
pzplot(zeroeshigh,poleshigh);



figure("Name","magnitude response of low pass filter");
h1 = DTFT(impulselow,w);
plot(w,abs(h1));
title("magnitude response of filter")
xlabel("frequency(radians)");
ylabel("magnitude");

figure("name", "impulse response of low pass filter")
stem(n,impulselow);

figure("Name","magnitude response of high pass filter");
h1 = DTFT(impulsehigh,w);
plot(w,abs(h1));
title("magnitude response of filter")
xlabel("frequency(radians)");
ylabel("magnitude");

figure("name", "impulse response of high pass filter")
stem(n,impulsehigh);

%% images
figure("Name","cat")
y = 0:1199;
x = 0:844;
img = imread('cat.jpg');
imshow(img);
%% pass through low pass filter
figure("Name","badcatlow");
badcat = convn(double(img), impulselow, 'same');
imshow(rescale(real(badcat)));
%% pass through high pass filter
figure("Name","badcathigh");
badcat = convn(double(img), impulsehigh, 'same');
imshow(rescale(real(badcat)));
%% end images
%% ALL FUNCTIONS SUPPORTING THIS CODE %%
function pzplot(b,a)
% PZPLOT(B,A)  plots the pole-zero plot for the filter described by
% vectors A and B.  The filter is a "Direct Form II Transposed"
% implementation of the standard difference equation:
% 
%    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 

    % MODIFY THE POLYNOMIALS TO FIND THE ROOTS 
    b1 = zeros(max(length(a),length(b)),1); % Need to add zeros to get the right roots
    a1 = zeros(max(length(a),length(b)),1); % Need to add zeros to get the right roots
    b1(1:length(b)) = b;    % New a with all values
    a1(1:length(a)) = a;    % New a with all values

    % FIND THE ROOTS OF EACH POLYNOMIAL AND PLOT THE LOCATIONS OF THE ROOTS
    h1 = plot(real(roots(a1)), imag(roots(a1)));
    hold on;
    h2 = plot(real(roots(b1)), imag(roots(b1)));
    hold off;

    % DRAW THE UNIT CIRCLE
    circle(0,0,1)
    
    % MAKE THE POLES AND ZEROS X's AND O's
    set(h1, 'LineStyle', 'none', 'Marker', 'x', 'MarkerFaceColor','none', 'linewidth', 1.5, 'markersize', 8); 
    set(h2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor','none', 'linewidth', 1.5, 'markersize', 8); 
    axis equal;
    
    % DRAW VERTICAL AND HORIZONTAL LINES
    xminmax = xlim();
    yminmax = ylim();
    line([xminmax(1) xminmax(2)],[0 0], 'linestyle', ':', 'linewidth', 0.5, 'color', [1 1 1]*.1)
    line([0 0],[yminmax(1) yminmax(2)], 'linestyle', ':', 'linewidth', 0.5, 'color', [1 1 1]*.1)
    
    % ADD LABELS AND TITLE
    xlabel('Real Part')
    ylabel('Imaginary Part')
    title('Pole-Zero Plot')
    
end


function circle(x,y,r)
% CIRCLE(X,Y,R)  draws a circle with horizontal center X, vertical center
% Y, and radius R. 
%
    
    % ANGLES TO DRAW
    ang=0:0.01:2*pi; 
    
    % DEFINE LOCATIONS OF CIRCLE
    xp=r*cos(ang);
    yp=r*sin(ang);
    
    % PLOT CIRCLE
    hold on;
    plot(x+xp,y+yp, ':', 'linewidth', 0.5, 'color', [1 1 1]*.1);
    hold off;
    
end



function H = DTFT(x,w)
% DTFT(X,W)  compute the Discrete-time Fourier Transform of signal X
% acroess frequencies defined by W. 

    H = zeros(length(w),1);
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.'*(nn-1));
    end
    
end


function xs = shift(x, s)
% SHIFT(x, s) shifts signal x by s such that the output can be defined by 
% xs[n] = x[n - s]

    % INITIALIZE THE OUTPUT
    xs = zeros(length(x), 1);
    
    for n = 1:length(x)
        % CHECK IF THE SHIFT IS OUT OF BOUNDS FOR THIS SAMPLE
        if n-s > 0 && n-s < length(x)
            % SHIFT DATA
            xs(n) = x(n-s);
        end
    end

end

function [b,a] = butter(N, Wn)
%  NOTE: JUST REPEATING MATLAB'S BUTTER INSTRUCTIONS HERE... 
%  THIS FUNCTION ONLY HAS A SIMPLE LOW PASS FUNCTIONALITY
%BUTTER Butterworth digital and analog filter design.
%   [B,A] = BUTTER(N,Wn) designs an Nth order lowpass digital
%   Butterworth filter and returns the filter coefficients in length
%   N+1 vectors B (numerator) and A (denominator). The coefficients
%   are listed in descending powers of z. The cutoff frequency
%   Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to
%   half the sample rate.
%
%   If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an
%   order 2N bandpass filter with passband  W1 < W < W2.
%   [B,A] = BUTTER(N,Wn,'high') designs a highpass filter.
%   [B,A] = BUTTER(N,Wn,'low') designs a lowpass filter.
%   [B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
%
%   When used with three left-hand arguments, as in
%   [Z,P,K] = BUTTER(...), the zeros and poles are returned in
%   length N column vectors Z and P, and the gain in scalar K.
%
%   When used with four left-hand arguments, as in
%   [A,B,C,D] = BUTTER(...), state-space matrices are returned.
%
%   BUTTER(N,Wn,'s'), BUTTER(N,Wn,'high','s') and BUTTER(N,Wn,'stop','s')
%   design analog Butterworth filters.  In this case, Wn is in [rad/s]
%   and it can be greater than 1.0.
%
%   % Example 1:
%   %   For data sampled at 1000 Hz, design a 9th-order highpass
%   %   Butterworth filter with cutoff frequency of 300Hz.
%
%   Wn = 300/500;                   % Normalized cutoff frequency
%   [z,p,k] = butter(9,Wn,'high');  % Butterworth filter
%   [sos] = zp2sos(z,p,k);          % Convert to SOS form
%   h = fvtool(sos);                % Plot magnitude response
%
%   % Example 2:
%   %   Design a 4th-order butterworth band-pass filter which passes
%   %   frequencies between 0.15 and 0.3.
%
%   [b,a]=butter(2,[.15,.3]);        % Bandpass digital filter design
%   h = fvtool(b,a);                 % Visualize filter
%
%   See also BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ,
%   FILTER, DESIGNFILT.

%   Author(s): J.N. Little, 1-14-87
%          J.N. Little, 1-14-88, revised
%          L. Shure, 4-29-88, revised
%          T. Krauss, 3-24-93, revised

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.


% % Poles
% ptemp = exp(1i*(pi*(1:2:n-1)./(2*n) + pi/2));
% 
% % ASSIGN POLES
% p = zeros(2*len+1,1); p(end) = -1;
% for k = 1:len
%     p(2*k-1) = ptemp(k);
%     p(2*k) = conj(ptemp(k));
% end

W0 = 0;
fs = 1;
Wn = Wn*pi;
%N = 4;

% COMPUTE Z-TRANSFORM POLES FROM LAPLACE POLES
% sp = zeros(N,1);
% zp = zeros(N,1);
% for k = 1:floor(N/2)
%     sp(k) = 1j*W0+Wn*exp(1j*(pi*(2*k-1)/(N) + pi/2));
%     zp(k) = (1+(sp(k)/2/fs))./(1-(sp(k)/2/fs));
% end
% for k = 1:ceil(N/2)
%     sp(k+floor(N/2)) = -1j*W0+Wn*exp(1j*(pi*(2*k-1)/(N) + pi/2));
%     zp(k+floor(N/2)) = (1+(sp(k+floor(N/2))/2/fs))./(1-(sp(k+floor(N/2))/2/fs));
% end

sp = zeros(N,1);
for k = 1:N
    sp(k) = Wn*exp(1j*(pi*(2*k-1)/(2*N) + pi/2));
end
[bs,as] = pz2ba(sp,zeros(length(sp),1));
[r,p,r0] = residue(bs,as);

r0 = 0;

% COMPUTE Z-TRANSFORM POLE LOCATIONS
zp = zeros(N,1);
for k = 1:N
    zp(k) = exp( p(k) );
end

% RECOMBINE
mypoly = zeros(N+1,1);
for kk = 1:N
    tmpkk = 1;
    for rr = 1:N
        if rr ~= kk
            tmpkk = conv(tmpkk, [1 zp(rr)]);
        end
    end
    mypoly = mypoly + [0; r(kk)*tmpkk(:)];
end
tmpkk = 1;
for rr = 1:N
    tmpkk = conv(tmpkk, [1 zp(rr)]);
end
mypoly = mypoly + r0*tmpkk(:);

mypoly = real(mypoly);

mydenom = 1;
for rr = 1:N
    mydenom = conv(mydenom, [1 -zp(rr)]);
end
mydenom = real(mydenom);

zz = roots(mypoly);
pp = roots(mydenom);

%[bs,as] = ba2pz(mypoly,zeros(length(sp),1));



% GET B AND A COEFFICIENTS
[b,a] = pz2ba(pp,zz);
b = b.*prod(1-pp)./prod(1-zz);         % Normalize Gain
b = real(b);
a = real(a);

end




function [b,a] = pz2ba(p,z)
% PZ2BA(P,Z)  Converts poles P and zeros Z to filter coefficients
%             B and A
%
% Filter coefficients are defined by:
%    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 

    % CONVERT ROOTS (POLES AND ZEROS) INTO POLYNOMIALS
    b = poly(z);
    a = poly(p);

end

function [p,z] = ba2pz(b,a)
% BA2PZ(B,A)  Converts filter coefficients B and A into poles P and zeros Z
% 
% Filter coefficients are defined by:
%    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 

    % MODIFY THE POLYNOMIALS TO FIND THE ROOTS 
    b1 = zeros(max(length(a),length(b)),1); % Need to add zeros to get the right roots
    a1 = zeros(max(length(a),length(b)),1); % Need to add zeros to get the right roots
    b1(1:length(b)) = b;    % New a with all values
    a1(1:length(a)) = a;    % New a with all values

    % FIND THE ROOTS OF EACH POLYNOMIAL
    p = real(roots(a1))+1j*imag(roots(a1));
    z = real(roots(b1))+1j*imag(roots(b1));
    
end


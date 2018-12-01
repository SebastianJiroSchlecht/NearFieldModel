function [num,den,fs] = nearFieldModel(theta, rho_n)
%nearFieldModel - computes the IIR filter for the near field transfer
%function as described in S. Spagnol, E. Tavazzi, and F. Avanzini,
%"Distance rendering and perception of nearby virtual sound sources with a 
%near-field filter model,"
%in Applied Acoustics, vol. 115, pp. 61-73, Jan. 2017.
%
% Computes the 
%
% Syntax:  [num,den,fs] = nearFieldModel(theta, rho_n)
%
% Inputs:
%    theta - Angle of incidence in radiants
%    rho_n - Distance of source, normalized to head radius
%
% Outputs:
%    num - Numerator of model filter
%    den - Denominator of model filter
%    fs  - Sampling frequency used
% Example: 
%    [num,den,fs] = nearFieldModel(0.3, 2.1)
%
% Other m-files required: clip, speedOfSound
%
% See also: visualize_nearFieldModel.m
% Author: Dr.-Ing. Sebastian Jiro Schlecht, 
% International Audio Laboratories, University of Erlangen-Nuremberg
% email address: sebastian.schlecht@audiolabs-erlangen.de
% Website: sebastianjiroschlecht.com
% 21. November 2018; Last revision: 21. November 2018

% constants
fs = 48000;
a_0 = 8.75 / 100;
a = 8.75 / 100;


% Table 1
data = [0  12.97  -9.69  -1.14 0.219  -4.39 2.123  -0.55    -0.06 0.457 -0.67 0.174 -1.75 0.699;
    10  13.19 234.2 18.48  -8.5  -4.31  -2.78 0.59          -0.17 0.455 0.142 -0.11 -0.01 -0.35;
    20  12.13  -11.2  -1.25 0.346  -4.18 4.224  -1.01       -0.02 -0.87 3404  -1699 7354  -5350;
    30  11.19  -9.03  -1.02 0.336  -4.01 3.039  -0.56       -0.32 0.465 -0.91 0.437 -2.18 1.188;
    40  9.91  -7.87  -0.83 0.379  -3.87  -0.57 0.665        -1.13 0.494 -0.67 0.658 -1.2  0.256;
    50  8.328  -7.42  -0.67 0.421  -4.1  -34.7 11.39        -8.3  0.549 -1.21 2.02  -1.59 0.816;
    60  6.493  -7.31  -0.5 0.423  -3.87 3.271  -1.57        0.637 0.663 -1.76 6.815 -1.23 1.166;
    70  4.455  -7.28  -0.32 0.382  -5.02 0.023  -0.87       0.325 0.691 4.655 0.614 -0.89 0.76;
    80  2.274  -7.29  -0.11 0.314  -6.72  -8.96 0.37        -0.08 3.507 55.09 589.3 29.23 59.51;
    90  0.018  -7.48  -0.13 0.24  -8.69  -58.4 5.446        -1.19 -27.4 10336 16818 1945  1707;
    100   -2.24  -8.04 0.395 0.177  -11.2 11.47  -1.13      0.103 6.371 1.735 -9.39 -0.06 -1.12;
    110   -4.43  -9.23 0.699 0.132  -12.1 8.716  -0.63      -0.12 7.032 40.88 -44.1 5.635 -6.18;
    120   -6.49  -11.6 1.084 0.113  -11.1 21.8  -2.01       0.098 7.092 23.86 -23.6 3.308 -3.39;
    130   -8.34  -17.4 1.757 0.142  -11.1 1.91 0.15         -0.4  7.463 102.8 -92.3 13.88 -12.7;
    140   -9.93  -48.4 4.764 0.462  -9.72  -0.04 0.243      -0.41 7.453 -6.14 -1.81 -0.88 -0.19;
    150   -11.3 9.149  -0.64  -0.14  -8.42  -0.66 0.147     -0.34 8.101 -18.1 10.54 -2.23 1.295;
    160   -12.2 1.905 0.109  -0.08  -7.44 0.395  -0.18      -0.18 8.702 -9.05 0.532 -0.96 -0.02;
    170   -12.8  -0.75 0.386  -0.06  -6.78 2.662  -0.67     0.05  8.925 -9.03 0.285 -0.9  -0.08;
    180   -13  -1.32 0.45  -0.05  -6.58 3.387  -0.84        0.131 9.317 -6.89 -2.08 -0.57 -0.4];

% assign variables
alpha = data(:,1);
p11 = data(:,2);
p21 = data(:,3);
q11 = data(:,4);
q21 = data(:,5);
p12 = data(:,6);
p22 = data(:,7);
q12 = data(:,8);
q22 = data(:,9);
p13 = data(:,10);
p23 = data(:,11);
p33 = data(:,12);
q13 = data(:,13);
q23 = data(:,14);

% convert and clip input arguments
theta = rad2deg(abs(theta));
rho = clip(rho_n, [1.25 Inf]);

% Eq (8), (13) and (14)
G_0 = ( (p11.*rho + p21) ./ (rho.^2 + q11.*rho + q21) ) ;
G_inf = ( (p12.*rho + p22) ./ (rho.^2 + q12.*rho + q22) ) ;
f_c = (p13.*rho.^2 + p23.*rho + p33) ./ (rho.^2 + q13.*rho + q23);

% denormalize
f_c = f_c*speedOfSound / (2*pi*a);

% linear interpolate at theta
iG_0 = interp1(alpha, G_0, theta);
iG_inf = interp1(alpha, G_inf, theta);
if_c = interp1(alpha, f_c, theta);

% Eq. (12), (10), and (11)
V_0 = db2mag(iG_inf);
tanF = tan(pi*a_0/a*if_c/fs);
a_c = (V_0 .* tanF - 1)./(V_0 .* tanF + 1);

% Eq (10)
V = (V_0 - 1) / 2;
b0 = V.*(1 - a_c) + 1;
b1 = V.*(a_c - 1) + a_c;
a1 = a_c;
a0 = ones(size(a1));

num = db2mag(iG_0) *  [b0 b1];
den = [a0, a1];

global m nfev Fvals Xhist Fhist Deltahist delta % For unedf_f checkpointing
addpath('./src/poundersb');
addpath('./src/poundersb/minq5/');

% SLy4 to UNEDF0 run bounds
Bl = [.15, -16.2, 190, 28, 40, .9, -Inf(1,6)];
Bu = [.17, -15.8, 230, 36, 100, 1.5, Inf(1,6)];
% Scaling bounds:
L = [.14 -17 170 27 30 .8 -200 -200 -300 -300 -200 -200];
U = [.18 -15 270 37 70  2  100  100 -100 -100  150  150];

Unedf_0= [
0.1605260d0
-16.05590d0
 230.00d0
 30.54290d0
 45.08040d0
 0.90d0
 -55.26060d0
-55.62260d0
-170.3740d0
-199.2020d0
-79.53080d0
 45.63020d0
]';

Unedf_1ex = [
 0.158366731281203
 -15.8000000000000
  220.000000000000
  28.3839520550152
  40.0000000000000
  1.00187172237670
 -44.6016362024908
 -180.956465344078
 -187.468585444369
 -207.209421806994
 -74.3391307279475
 -38.8371791097226
   0.813550803208000
 ]'; %'

Unedf_1 = [0.158706769332587
-1.580000000000000e+01
2.200000000000000e+02
28.9867890577721
40.0047904804136
0.992423332283364
-45.1351310222373
-145.382167908057
-186.065399575124
-206.579593890106
-74.0263331764599
-35.6582611147917]';

SLy4 = [0.159538756711733
-15.972149141444596
229.900964482603626
32.004302815052007
45.961751480461579
1.439546988976078
-76.996203124999993
15.657135124999999
-285.84
-285.84
-92.250000000000000
-30.750000000000000
]'; %'

DME_LO_new = [
0.1640175211422964
-15.97328699890108
190.000
32.67940306677653
40.25000000000000
0.9017357885507260
-44.90707953241770
19.61223611419416
-181.8690019064476
-210.6411340785048
-85.43134376981442
-80.28050191487668 ]'; %'

DME_NLO_new = [
0.1576434284368679
-15.8000000000000
230.0000000000000
30.14869224442319
40.00000000000000
1.140201102344124
-38.76206161053054
-97.77475488543688
-213.1893875314728
-227.6904159491677
-81.47618598693505
-40.11569763601180 ]'; %'

DME_N2LO = [
0.155053
-15.9816
219.2050
33.8459
45.8224
1.2295
-190.1675
58.2696
-227.5362
-229.4380
-184.0530
-34.7645 ]'; %'

func = @unedf_f; % [f h] Handle so func(x) returns 1-by-m vector F
n = 12;         % [int] Dimension (number of continuous variables)
npmax = .5*(n+1)*(n+2);  % [int] Max # of interp. points (>n+1) (2*n+1)
nfmax = 100;     % [int] Max # of fun. evals done by pounders (>n+1) (100)
printf = true;  % [log] 1 Indicates you want output to screen (1)
gtol = 1e-10;   % [dbl] Tolerance for 2-norm of model gradient (>0) (1e-4)
delta = .1;     % [dbl] Positive trust region radius (.1)
X0 = (DME_LO_new-L)./(U-L);      % [dbl] [max(nfs,1)-by-n] Set of initial pts (zeros(1,n))
xind = 1;       % [int] Index of point in X0 at which to start from (1)
nfs = 0;        % [int] Number of f values (at X0) known in advance (0)
m = 115;        % [int] Number of residuals
F0 = [];        % [dbl] [nfs-by-m] Set of known function values  ([])
Low = (Bl-L)./(U-L); % [dbl] [1-by-n] Vector of lower bounds (-Inf(1,n))
Upp = (Bu-L)./(U-L); % [dbl] [1-by-n] Vector of upper bounds (Inf(1,n))

% Choose your subproblem solver:
global spsolver
spsolver=2;  %Arnold Neumaiers minq

%X0=X0(1:12);
%Low = Low(1:12);
%Upp = Upp(1:12);
npmax = (n+1)*(n+2)/2;

% Checkpointing initialization:
Fvals = zeros(nfmax,1); Xhist = zeros(nfmax,n); Fhist = zeros(nfmax,m);
Deltahist = zeros(nfmax,1); nfev = 0;

[X,F,flag,xkin] = ...
    pounders(func,X0,n,npmax,nfmax,gtol,delta,nfs,m,F0,xind,Low,Upp,printf);
save(['Poundersb_test'],'X','F','flag','xkin')

function Fvec = unedf_f(X)
% This routine is used to specify a generic function called by POUNDerS
% Takes n-by-1 column vector X and outputs m-by-1 vector Fvec.


% Global vars used in pounders_driver
global m nfev Fvals Xhist Fhist Deltahist delta
n = length(X); % Problem dimension
nfev = nfev +1;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Generate the vector through internal or external routine between the ++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Remove existing hel files
if 1==0
    unix('rm -f *.hel');
end


if 1==0 % Use some other hel files
    unix('cp helrepo/*.hel .');
end

% Write x values to files
L = [.14 -17 170 27 30 .8 -200 -200 -300 -300 -200 -200];%-100 -100];
U = [.18 -15 270 37 70  2  100  100 -100 -100  150  150];%100 100];
fid = fopen('x.in','wt');

fprintf(fid,'%1.15e\n',[X.*(U-L)+L,1]); %X.*(U-L)+L
fclose(fid);

% Run the HFBTHO script
unix('./runc.sh');

% Read in f from file:
fid = fopen('fval.out','r');
Fvec = fscanf(fid,'%e');
fclose(fid);

% Clean up:
%delete('x.in'); delete('fval.out');

Fvec=Fvec'; % Make into row vector

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% [OPTIONAL] CORRECTION DONE TO AVOID HUGE VALUES:
Fvec = max(Fvec,-1e64);
Fvec = min(Fvec,1e64);

% Save details for the final output:
normf = sum(Fvec.^2)
fprintf('||f||^2 = %f\n',normf);
Fvals(nfev) = normf; % Stores the history of the objective vals
Fhist(nfev,:) = Fvec;      % Stores the history of the fvec vals
Xhist(nfev,:) = X;         % Stores the history of the x vals
Deltahist(nfev) = delta;
save(strcat('H_checkpoint_',date),'Fvals','Xhist','Fhist','Deltahist');

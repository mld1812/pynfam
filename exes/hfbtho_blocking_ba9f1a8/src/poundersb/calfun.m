function F=calfun(x)

x=x(:); % Turn into column vector
n=length(x);
%rand('state',0);
%A=rand(n);
%a=1; % Can change this to be many things
%y=norm(A*(x-a*ones(n,1)))^2 +norm(x)^4;
%y=sqrt(y);


F = [x', .01*ones(1,108-n)]; %x';%+(x.^2)';
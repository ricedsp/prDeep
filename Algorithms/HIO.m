%Basic implementation of HIO algorithm
%Fienup, James R. "Phase retrieval algorithms: a comparison." Applied optics 21.15 (1982): 2758-2769.
%Adapted from work by Prasanna Rangarajan
function x = HIO( y, A, At, support, beta, iters, x_init )
%Assumes x is real
if isnumeric(A)
    A=@(x) A*x(:);
    At=@(z) A'*z(:);
end

% random initial guess
if nargin < 7
  x = rand(size(At(y))); % some random real non-negative signal
  x(1)=1;
else
  x = x_init;
end

Pm=@(x) At(y.*exp(1i*angle(A(x))));

%Main loop
for k = 1:iters
	Pmx=Pm(x);
    inds=logical(isreal(Pmx).*support.*(real(Pmx)>0));
	x(inds)=Pmx(inds);
    x(~inds)=x(~inds)-beta*Pmx(~inds);
end
end  
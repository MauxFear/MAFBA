function model = setWeights_rand(model,group,opt,multiplier,shape)
% function model = setWeights_rand(model,group,opt,multiplier,shape)
%
% 
if ~isfield(model,'w') || ~isfield(model,'protGroup')
    error('Error: no ''w'', or ''protGroup'' field in the model.');
end;

if ~exist('group','var')
    group=0; % all rxns
elseif group<=0 || group>length(model.protGroup)
    error('group>length(model.protGroup) in setWeights_rand()');
end

if ~exist('shape','var') % (delta)
	disp('Error: ''shape'' parameter missing. Exiting.');
    return;
end

%%%%%%%%%%%% Extraction of the numbers

rr=model.protGroup(group).rxns; 
Nrr=length(rr);

if     strcmpi(opt,'gamma')
    % Option 1: extract the weights from a gamma distribution
    model.w(rr)=multiplier*generateGammaRnd(Nrr,shape);
    
elseif strcmpi(opt,'expbox')
    % Option 2: extract the weights from a x=10^unif. pdf, that is, p(x) ~ 1/x
    model.w(rr)=multiplier*generateExpBoxRand(Nrr,shape);
    
elseif strcmpi(opt,'lognormal')
    % Option 3: extract the weights from a lognormal distribution
    model.w(rr)=multiplier*generateLogNormalRand(Nrr,shape);
else
    error('Unrecognized ''opt'' variable in setWeights_rand().');
end

end

function x=generateGammaRnd(N,k)
%generateGammaRnd(N,k)
% Generates N random numbers according to a Gamma distribution
%  f(x|k,theta) ~ x^(k-1) e^(-x/theta) ~ gampdf(x,k,theta)
%
% The shape parameter k>0 is chosen by the user, while
% the scale parameter b is set such that the expected value
% of the random variabile is 1 for any a:
%
% E(x)=1 
% V(x)= 1/k
%
% (tends to a gaussian for k->\infty)
%
if ~exist('k','var') || k <= 0
    k=1; % default value
end
theta=gamma(k)/gamma(k+1);

x= gamrnd(k,theta,[N,1]);
end

function x=generateExpBoxRand(N,w)
% Generates a random variable distributed as the exponential of an
% uniform distribuited number. If p(x)=1/w for x in [a,a+w], and calling
% y=c^x, we get the distribution:
% p(y)=1/(ln(c)*w*y) in the interval [c^a,c^(a+w)].
% The expected value is:
%
% <y>=c^a/(ln(c)*w)*(c^w-1)
%
% Choosing <y>=1, we find a in terms of c and w:
%
% a=ln(ln(c)*w/(c^w-1))/ln(c)
%
if w==0
    x=ones(N,1);
else
    c=10; % If c=10, w is the number of decades spanned by y (delta)
    a=log(log(c)*w/(c^w-1))/log(c);
    x= c.^(a + w*rand(N,1));
end
end

function x=generateLogNormalRand(N,sigma)
% Extracts N numbers from a lognormal distribution

if sigma==0
    x=ones(N,1);
else
    x=zeros(N,1);
    for i=1:N
        x(i)=10^( normrnd(0,sigma));
    end
end
end





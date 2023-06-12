%% Code to enumerate the exact results derived in "Exact sharp-fronted solutions
%% for nonlinear diffusion on evolving domains" by Johnston and Simpson, 2023.
%% The default parameters here generate the red solution profiles in Figure 
%% 1(a)-(b) in the main manuscript.

close all
clear all

%% Define model parameters.

L0 = 50;                                    % Define initial domain size.
beta = 0.01;                                % Define domain growth rate.
m = 1.5;                                    % Define order of diffusivity (m*c(x,t)^{m-1}).
n = 1;                                      % Define dimensionality.
c = 1;                                      % Define free parameter related to initial mass.
domainGrowthType = 'linear';                % Choose growth type 'linear' or 'exponential'.
alpha = n/(n*(m-1)+2);                      % Parameter in porous medium equation solution.
gamma = 1/(n*(m-1)+2);                      % Parameter in porous medium equation solution.
tEnd = 500;                                 % Final time for the solution profile to be evaluated.

%% Define solution components

if strcmpi(domainGrowthType,'linear')
    domainSize = @(t) L0 + beta*t;          % Domain growth function
    if beta > 0
        timeTransform = @(t) n^(m+1)*(L0 - L0*(L0./(L0+beta*t)).^m)/(beta*m);   % Time transformation
    elseif beta == 0
        timeTransform = @(t) n^(m+1)*t;                                         % Time transformation
    else
        error('Positive growth rates only')                                     % Only supports positive growth rates
    end  
    % Point where support intersects with moving boundary. If imaginary, does not occur.
    boundaryCollision = L0 * (1 - beta*m*n^(-m-1)/L0*(gamma*(m-1)*L0^2/(2*c*m))^(1/(2*gamma)))^(-1/m);  
    betaCritical = (m*n^(-m-1)/L0*(gamma*(m-1)*L0^2/(2*c*m))^(1/(2*gamma)))^-1;                         % Critical value of beta
    tStar = (boundaryCollision - L0)/beta;                                                              % Critical value of time
elseif strcmpi(domainGrowthType,'exponential')
    domainSize = @(t) L0*exp(beta*t);
    if beta > 0
        timeTransform = @(t) n^(m+1)*(1-exp(-beta*(m+1)*t))/(beta*(m+1));       % Time transformation
    elseif beta == 0
        timeTransform = @(t) n^(m+1)*t;                                         % Time transformation
    else
        error('Positive growth rates only')                                     % Only supports positive growth rates
    end
    % Point where support intersects with moving boundary. If imaginary, does not occur.
    boundaryCollision = L0 * (1 - beta*(m+1)*n^(-m-1)*(gamma*(m-1)*L0^2/(2*c*m))^(1/(2*gamma)))^(-1/(m+1));  
    betaCritical = ((m+1)*n^(-m-1)*(gamma*(m-1)*L0^2/(2*c*m))^(1/(2*gamma)))^-1;                            % Critical value of beta
    tStar = log(boundaryCollision/L0)/beta;                                                                  % Critical value of time
end

domainRatio = @(t) L0./domainSize(t);                                       % Ratio of initial to current domain size.
nTimePoints = 1001;                                                         % Number of time points at which to plot results.
time = [1e-4,linspace(tEnd/nTimePoints,tEnd,nTimePoints-1)];                % Time values at which to plot results.
nDomainPoints = 1001;                                                       % Number of spatial points at which to plot results.

solution = zeros(nTimePoints,nDomainPoints);                                % Initialise solution profile storage.
domainPoints = zeros(nTimePoints,nDomainPoints);                            % Initialise domain points storage.
zeroBoundaryLocation = zeros(nTimePoints,1);                                % Support width (based on solution profile).

for iTime = 1:nTimePoints
    domainBoundary = domainSize(time(iTime));                                                                           % Current location of domain boundary.
    domainPoints(iTime,:) = linspace(-domainBoundary,domainBoundary,nDomainPoints);                                     % Current points at which to evaluate solution.
    tmpSolutionComponent = c-(m-1)/(2*m)*gamma*(domainPoints(iTime,:)*domainRatio(time(iTime))).^2 ...
        * timeTransform(time(iTime))^(-2*gamma);                                                                        % Part of PME solution.
    tmpSolutionComponent(tmpSolutionComponent<0) = 0;
    solution(iTime,:) = n*timeTransform(time(iTime))^(-alpha)*tmpSolutionComponent.^(1/(m-1))*domainRatio(time(iTime)); % Transformed PME solution.
    zeroLocations = find(solution(iTime,:)==0);
    zeroBoundary = min(zeroLocations(zeroLocations>ceil(nDomainPoints/2)));
    % Find width of support from solution profile.
    if isempty(zeroBoundary)
        zeroBoundary = nDomainPoints;
    end
    zeroBoundaryLocation(iTime) = domainPoints(iTime,zeroBoundary);
end

supportWidth = sqrt(2*c*m/((m-1)*gamma)*timeTransform(time).^(2*gamma))./domainRatio(time); % Analytic expression of support width over time

% Plot support width over time (option to check against width obtained from
% solution profile.
figure(1); plot(time,supportWidth,'linewidth',3); hold on;
%hold on; plot(time,zeroBoundaryLocation)

% Plot final solution profile.
figure(2); plot(domainPoints(end,:),solution(end,:),'linewidth',3); hold on;
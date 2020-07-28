% Le Xuan Thang 
% Cuckoo search (CS)
% 26-06-2019
% University of Transport and Communications
% Project: Modal updating problem
%------------------------------------------------------------%


%% empty space 
clc
clear

%% start
tic;

%% Function
CostFunction=@(x)MyThuan34m(x);

%% Parameters
MaxIt=100; % Maximum iteration
nPop=10;   % Population size
Boundary=[2.96E10 4.1884E10;
            1.85E11 2.3E11
            2.856E10 3.687E10
            2.874E10 4.064E10
            2.944E10 4.01E10
            2.874E10 4.01E10
            2.944E10 4.01E10
            2.874E10 3.71E10
            2.821E10 4.01E10];
nVar=length(Boundary);  % Number of variables
VarSize=[1 nVar];       % Size of variables
VarMin=Boundary(:,1);   % Lower boundary
VarMax=Boundary(:,2);   % Upper boundary
% Discovery rate of alien eggs/solutions
pa=0.25;

%% Structure
GlobalBest.Cost = Inf;
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
particle=repmat(empty_particle,nPop,1);
for i = 1 : nPop
particle(i).Best.Cost=10^10;
particle(i).Best.Position=[];
end

%% Initialization
for i = 1 : nPop
   for j = 1 : nVar
    % Position
    particle(i).Position(1,j) = unifrnd(VarMin(j,1),VarMax(j,1));
   end
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position);
    
    % Update Particle Best Solution
    if particle(i).Cost < particle(i).Best.Cost
        particle(i).Best.Cost=particle(i).Cost;
        particle(i).Best.Position=particle(i).Position;
        if  particle(i).Best.Cost<GlobalBest.Cost
            GlobalBest=particle(i).Best;
        end
    end
end

%% Main Loop
it = 0;
BestCosts=zeros(MaxIt,1);
while it<MaxIt
    it=it+1;
    % Get new solution but keep current best
    % Levy exponent and coefficient
    % For details, see equation (2.21), Page 16 (chapter 2) of the book
    % X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/...
        (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    nest = vertcat (particle(:).Position);%vertcat connect two matrix.
    for i = 1 : nPop
        s=nest(i,:);
        u =randn(size(s))*sigma; %randn Normally distributed random numbers
        v=randn(size(s));
        step=u./abs(v).^(1/beta);%v is the same t in equation calculating CS; step is levy
        stepsize=0.01*step.*(s-GlobalBest.Position);
        % Now the actual random walks or flights
        s=s+stepsize.*randn(size(s)); 
        
        % Apply limit Position
        for j = 1:nVar
        s(1,j) = max(s(1,j),VarMin(j,1));
        s(1,j) = min(s(1,j),VarMax(j,1));
        end
        particle(i).New.Position=s;
        particle(i).New.Cost=CostFunction(particle(i).New.Position);
        if particle(i).New.Cost < particle(i).Best.Cost
            particle(i).Best=particle(i).New;
            if particle(i).Best.Cost<GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
    BestCosts(it,1)=GlobalBest.Cost;
    disp(['Iteration = ' num2str(it) '********** BestCost = ' num2str(BestCosts(it,1))])
    % Update the counter
      it=it+1; 
    % Discovery and randomization
    for i = 1:nPop
    nest(i,:) = vertcat(particle(i).Best.Position);
    end
    % Discovered or not -- a status vector
    K=rand(size(nest))>pa;
    % New solution by biased/selective random walks
    stepsize=rand*(nest(randperm(nPop),:)-nest(randperm(nPop),:));
    new_nest=nest+stepsize.*K;
    for i = 1:nPop
        s=new_nest(i,:);
        % Apply Limits
        for j = 1:nVar
        s(1,j) = max(s(1,j),VarMin(j,1));
        s(1,j) = min(s(1,j),VarMax(j,1));
        end
        particle(i).New.Position=s;
        % Evaluate
        particle(i).New.Cost=CostFunction(particle(i).New.Position);
        if particle(i).New.Cost < particle(i).Best.Cost
            particle(i).Best=particle(i).New;
            if particle(i).Best.Cost<GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
    BestCosts(it,1)=GlobalBest.Cost;
    disp(['Iteration = ' num2str(it) '********** BestCost = ' num2str(BestCosts(it,1))])
end

disp(['Vi tri tot nhat la' num2str(GlobalBest.Position) 'Gia tri tot nhat la' num2str(GlobalBest.Cost)])
BestSol=GlobalBest;


%% Results

toc;

% Plot
figure;
semilogy(BestCosts,'-b+','LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;



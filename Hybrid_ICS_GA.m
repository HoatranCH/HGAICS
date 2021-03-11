
clc;clear;
CostFunction=@(x) MyThuan34m(x);
tic;
%Boundary (apply your problem)
nVar=size(Boundary,1);

VarSize=[1 nVar];

VarMin=Boundary(:,1);
VarMax=Boundary(:,2);
MaxIt=100;
nPop=10;

%% CS Parameters
% Discovery rate of alien eggs/solutions
% pa=0.25;
pa_max=0.5;
pa_min=0.005;
%% GA parameters
pc=0.7;                 % Crossover Percentage, Phan tram lai cheo
nc=2*round(pc*nPop/2);  % Number of Offsprings (also Parnets), So phan tu duoc phep lai cheo
gamma_ga=0.4;              % Extra Range Factor for Crossover

pm=0.3;                 % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants
mu=0.1;         % Mutation Rate
beta=8;

%% Template of Population
% Population
empty_individual.Position=[];
empty_individual.Cost=[];

pop=repmat(empty_individual,nPop,1);

% Paticle
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

% Global Best
GlobalBest.Cost=inf;


% Khoi tao ngau nhien vi tri cua hat, vector,...
for i=1:nPop
    
    for j=1:nVar
    % Tao ngau nghien vi tri tu VarMin den VarMax cho VarSize phan tu
    pop(i).Position(1,j)=unifrnd(VarMin(j,1),VarMax(j,1)); 
    end
    % Tinh gia tri cua vi tri do
    pop(i).Cost=CostFunction(pop(i).Position(1,:));       
    
    % Nhap vi tri cho hat
    particle(i).Position=pop(i).Position;          
    
    % Tao ma tran van toc choi VarSize kich thuoc
    
    particle(i).Velocity=zeros(VarSize);            
    
    % Lay gia tri tinh duoc cua Ham cho tung gia tri cua hat
    particle(i).Cost=pop(i).Cost;                   
    
    % Cap nhat vi tri tot nhat cua Hat
    particle(i).Best.Position=particle(i).Position; % Vi tri tot nhat cua Hat
    particle(i).Best.Cost=particle(i).Cost;         % Gia tri tot nhat tuong ung cua hat
    
    % Cap nhat vi tri tot nhat toan cau
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;        % So sanh gia tri tot nhat cua Hat voi Nhom
        
    end
end




generation=0;    
tolerance=1;
while generation<MaxIt && tolerance>10^-12
generation=generation+1;
 %%   
Costs=[pop.Cost];                   % Dinh nghia gia tri

% Lay gia tri sap xep tu nho den lon theo cot, va thu tu sap xep cua no
[Costs, SortOrder]=sort(Costs);     
pop=pop(SortOrder);

% Lay gia tri tot nhat sau khi sap xep
BestSol=pop(1);                     

% Gia tri xau nhat
WorstCost=pop(end).Cost;        
    % Su dung CS tao ra dan so cua GA
    % Get new solution but keep current best
    % Levy exponent and coefficient
    % For details, see equation (2.21), Page 16 (chapter 2) of the book
    % X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
    beta_CS=3/2;
    sigma_CS=(gamma(1+beta_CS)*sin(pi*beta_CS/2)/(gamma((1+beta_CS)/2)*beta_CS*2^((beta_CS-1)/2)))^(1/beta_CS); %equation 15 in the paper
%     sigma_CS=0.6966;
    for i = 1:nPop
    nest(i,:) = vertcat(particle(i).Best.Position);
    end
    for i = 1 : nPop
        s=nest(i,:);
        u=randn(size(s))*sigma_CS;
        v=randn(size(s));
        step=u./abs(v).^(1/beta_CS); %equation 14 in the paper
        stepsize=0.01*step.*(s-GlobalBest.Position);
        % Now the actual random walks or flights
        % Update Position
        s=s+stepsize.*randn(size(s)); %equation 21 in the paper; randn(size(s)represents random walk (levy lamda); stepsize represents anphal(n)
         %%
        % Line 136 to line 142 is wrong or different with the way
        % calculating s of the paper or Yang, in the paper, stepsize =1,
        % and stepsize (in the paper) is the step in Matlab code; it means
        % line 142 is not necessary and is corrected as folows: s=s+step*alpha
        % Apply limit Position
        % Apply limit Position
        for j = 1:nVar
        s(1,j) = max(s(1,j),VarMin(j,1));
        s(1,j) = min(s(1,j),VarMax(j,1));
        end
        
        particle(i).New.Position=s;
        particle(i).New.Cost=CostFunction(particle(i).New.Position);
        if particle(i).New.Cost < particle(i).Best.Cost
            particle(i).Best=particle(i).New;
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
%     BestCost(generation,1)=GlobalBest.Cost;       
%     disp(['Iteration ' num2str(generation) ':---------Best Cost = ' num2str(BestCost(generation,1))]);
%     generation=generation+1;
    % Caculate pa
      pa=pa_max-(((pa_max-pa_min)/MaxIt)*generation);
    % Discovery and randomization
    for i = 1:nPop
    nest(i,:) = vertcat(particle(i).Best.Position);
    end
    % Discovered or not -- a status vector
    K=rand(size(nest))>pa;
    % New solution by biased/selective random walks
    stepsize=rand*(nest(randperm(nPop),:)-nest(randperm(nPop),:));
    % Update Position
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
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
        % Cap nhat tro lai cho dan so
    for i=1:nPop
        if particle(i).Best.Cost < pop(i).Cost
            pop(i).Position = particle(i).Best.Position;
            pop(i).Cost = particle(i).Best.Cost;
            if pop(i).Cost<GlobalBest.Cost
                GlobalBest.Cost=pop(i).Cost;
                GlobalBest.Position=pop(i).Position;
            end
        end
    end

%-------------------------------End Cuckoo Search-------------------------------------------%    
%%
%-----------------------------------Start GA -----------------------------------------------%

    P=exp(-beta*Costs/WorstCost);   % p(0)=e^(-beta*Costs/WorstCost)
    P=P/sum(P);  
    
    % Khoi tao kich thuoc ma tran population crossover theo hang
    pop_crossover=repmat(empty_individual,nc/2,2);
    % Chon 1 nua trong so phan tu duoc lai cheo
    for k=1:nc/2 
        % Chon ngau nhien the he cha me de lai giong
        % Select Parents Indices
        i1=RouletteWheelSelection(P);   % Gen cua cha
        i2=RouletteWheelSelection(P);   % Gen cua me
        p1=pop(i1);                     % Lay Gen cua cha
        p2=pop(i2);                     % Lay Gen cua cha
        
        % Apply Crossover
        % Lai cheo 
        [pop_crossover(k,1).Position, pop_crossover(k,2).Position]=Crossover(p1.Position,p2.Position,gamma_ga,VarMin,VarMax,nVar);
        
        % Evaluate Offsprings
        pop_crossover(k,1).Cost=CostFunction(pop_crossover(k,1).Position);
        pop_crossover(k,2).Cost=CostFunction(pop_crossover(k,2).Position);        
    end    
    
    % Lay tat ca gia tri vua roi xep theo 1 hang doc
    pop_crossover=pop_crossover(:);

    % Dot Bien
    pop_mutation=repmat(empty_individual,nm,1);
    
    for k=1:nm
        
        % Lua chon cha me
        i=randi([1 nPop]);
        p=pop(i);
        
        % Apply Mutation
        pop_mutation(k).Position=Mutate(p.Position,mu,VarMin,VarMax,nVar);
        
        % Evaluate Mutant
        pop_mutation(k).Cost=CostFunction(pop_mutation(k).Position);
        
    end
    


    % Gop toan bo dan so vao
    pop=[pop
        pop_crossover
        pop_mutation];  %ok
    
    % Sap xep lai dan so
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Tim gia tri xau nhat
    WorstCost=max(WorstCost,pop(end).Cost);
    
    % Cat Giam
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);
    

    % Tap hop nhung gia tri tot nhat tim duoc
    for i=1:nPop
%         Neu gia tri cua tap hop nho hon gia tri cua hat
        if pop(i).Cost < particle(i).Best.Cost
%             Cap nhat vi tri va gia tri cua hat
            particle(i).Best.Position = pop(i).Position;
            particle(i).Best.Cost = pop(i).Cost;
            if  particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
    
% calculating tolerance    
    if generation>100
    tolerance=abs(BestCost((generation-100),1)-GlobalBest.Cost);
    end
    
    BestCost(generation,1)=GlobalBest.Cost;       
    disp(['Iteration ' num2str(generation) ':---------Best Cost = ' num2str(BestCost(generation,1))]);
            

end
toc;
figure;
semilogy(BestCost,'-g*','LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

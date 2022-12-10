A = [1,0.05,0.07;0.05,1.05,0.15;0,0.2,1.65];
B = [0.31,0.1,0.25;0.14,0.58,0.5;0,0.2,0.5]; 
Q = [10,3,1;3,5,4;1,4,9]; 
R = 1*eye(3);
env = myDiscreteEnv(A,B,Q,R);
rng(0)

K0 = place(A,B,[0.4,0.8,0.5]);
agent = LQRCustomAgent(Q,R,K0);
agent.Gamma = 1;
agent.EstimateNum = 45;
trainingOpts = rlTrainingOptions(...
    'MaxEpisodes',45, ...
    'MaxStepsPerEpisode',500, ...
    'Verbose',true, ...
    'Plots','training-progress');
trainingStats = train(agent,env,trainingOpts);
simOptions = rlSimulationOptions('MaxSteps',200);
experience = sim(env,agent,simOptions);
totalReward = sum(experience.Reward)
[Koptimal,P] = dlqr(A,B,Q,R); 
x0 = experience.Observation.obs1.getdatasamples(1);
Joptimal = -x0'*P*x0;
rewardError = totalReward - Joptimal


len = agent.KUpdate;
err = zeros(len,1);
updates = agent.KUpdate

ks = agent.K;
kmag = norm(ks);
    



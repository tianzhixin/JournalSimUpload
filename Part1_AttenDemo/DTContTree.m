%% Date: 2017.08.24
%----------------------Simulation Results-----------------------------
% Successful
% Results saved in Kv_r_50Regions02Delta60s20170824.mat
%% K `matrix optimization with constraints of K's structure, tree stability, and greshgorin theorem.
% linear objective, obj=norm(K,'fro');

% Equality linear constraints, inequality elementwise linear constraints
% constraints=[constraints,K(i,j)==0];

% secondordercone linear constraints: 
% constraints=[constraints,Q_rowsum(:)<=1-delta];



%----------------------Problem Description-----------------------------
% Discrete Time
% Tree Network







%----------------------Problem setup-----------------------------
% load system_info
% load Asys_Bsys_3Regions
% load AsysBsys50Regions_Tree170824.mat

% A(2,3) = A(1,1)*0.01;
% A(3,2) = A(1,1)*0.02;

% r=0.1;
function [Kv_r]=DTContTree(T,delta,no_region,A_adj,Asys,Bsys,L)
% T=60;
% delta=0.20;  %constraints=[constraints,Q_rowsum(:)<=1-delta];
A = Asys;


%----------------------Yalmip varibles-----------------------------

no_beta=length(nonzeros(A_adj));  %A_adj: Adjacency matrix of tree T_ori with 1 diagonal
K=sdpvar(no_beta,no_region);
Q=sdpvar(no_region,no_region);
% Q_ori=sdpvar(no_region,no_region);
Q_rowsum=sdpvar(no_region,1);
R=sdpvar(no_region,no_region);

% obj = 0;
K_test = ones(no_beta,no_region);
constraints=[];

%----------------------Sparsity constraints on K Demonstration-------
for j=1:no_region
    i=1;
    for Ai=1:no_region
        for Aj=1:no_region
            if A_adj(Ai,Aj)~=0
                if (Ai==j)||(Aj==j)
                elseif (Ai==Aj)&&(A_adj(Ai,j)~=0)
                else
                   K_test(i,j)=0;
%                     obj=obj + k*(K(i,j)).^2;
                end
                i=i+1;
            end
        end
    end
end


%----------------------Sparsity constraints on K---------------------
for j=1:no_region
    i=1;
    for Ai=1:no_region
        for Aj=1:no_region
            if A_adj(Ai,Aj)~=0
                if (Ai==j)||(Aj==j)
                elseif (Ai==Aj)&&(A_adj(Ai,j)~=0)
                else
                    constraints=[constraints,K(i,j)==0];
%                     obj=obj + k*(K(i,j)).^2;
                end
                i=i+1;
            end
        end
    end
end

%----------------------Tree Stability Constraints---------------------
constraints=[constraints,Q==eye(no_region)+T*L*(Asys+Bsys*K)*inv(L)];
for i=1:no_region
    Q_rowsum(i)=sum(abs(Q(i,:)));
end
constraints=[constraints,Q_rowsum(:)<=1-delta];

%------------Lyapunov Stability Constraints on Modeled system--------------


% constraints=[constraints,[eye(no_region) Q; Q' eye(no_region)-r*eye(no_region)]>=0];

% constraints=[constraints,[eye(no_region) Q_ori; Q_ori' eye(no_region)-R]>=0];
% constraints=[constraints,R>=r*eye(no_region)];
% constraints=[constraints,R<=0.95*eye(no_region)];

%--------------------------Objective and Solver----------------------------

obj=norm(K,'fro'); %linear
% objective=norm(K); %linear

options = sdpsettings('verbose',1,'solver','gurobi');  % gurobi: QP solver 
% options = sdpsettings('verbose',1,'solver','sedumi'); %sedumi: SDP solver

optimize(constraints,obj,options)

%% The values of the variables
% r

    Kv=value(K);
    Qv=value(Q);

    Qv_eig=eig(Qv);
    Unstable_eig_Q = find(abs(Qv_eig)>=1);
    if Unstable_eig_Q>0
        Unstable_eig_Q
    end
%     Kv_froNorm = norm(Kv,'fro');
%     Kv_norm=norm(Kv);
%     Kv_1norm=norm(Kv,1)

%% Constraints satisfied?
% for j=1:no_region
%     i=1;
%     for Ai=1:no_region
%         for Aj=1:no_region
%             if A_adj(Ai,Aj)~=0
%                 if (Ai==j)||(Aj==j)
%                 elseif (Ai==Aj)&&(A_adj(Ai,j)~=0)
%                 else
%                     if Kv(i,j)>=10e-8
%                         error = 1
%                         i
%                         j
%                     end
%                     
%                 end
%                 i=i+1;
%             end
%         end
%     end
% end
% Qv-eye(no_region)-T*(Asys+Bsys*Kv)
% Qv_rowsum=value(Q_rowsum);
% all(Qv_rowsum<1)
% all(eig([eye(no_region) Qv; Qv' eye(no_region)-r*eye(no_region)]) > 0)
%% Kv implimentation
% Although constraints=[constraints,K(i,j)==0]; are enforced, Kv does not
% have strict zero elements. However, they are small.
% To implement Kv in simulink models, Kv_r is to round Kv near zero
% elements to zero.

ZeroEle = (abs(Kv)<10e-10);
SparMat = -1*(ZeroEle-1);
Kv_r = Kv.*SparMat;

% Test whether Kv_r satisfies the equality elementwise constraints.
for j=1:no_region
    i=1;
    for Ai=1:no_region
        for Aj=1:no_region
            if A_adj(Ai,Aj)~=0
                if (Ai==j)||(Aj==j)
                elseif (Ai==Aj)&&(A_adj(Ai,j)~=0)
                else
                    if abs(Kv_r(i,j))>=10e-10
                        error = 1
                        i
                        j
                    end
                    
                end
                i=i+1;
            end
        end
    end
end
end

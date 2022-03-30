%% Composed Matrices

[Aol12,BF12,Be12,Bf12,Ceta12,Detae12,Detaf12,Cz12,Dze12,Dzf12,~,~] = getN(SS.A(:,:,1),SS.Bu(:,:,1),SS.Bw(:,:,1),SS.Cy(:,:,1),SS.Dyv(:,:,1),eigL(2),a,b,true,false,wnoise);
[Aol1N,BF1N,Be1N,Bf1N,Ceta1N,Detae1N,Detaf1N,Cz1N,Dze1N,Dzf1N,~,~] = getN(SS.A(:,:,1),SS.Bu(:,:,1),SS.Bw(:,:,1),SS.Cy(:,:,1),SS.Dyv(:,:,1),eigL(N),a,b,true,false,wnoise);
[Aol22,BF22,Be22,Bf22,Ceta22,Detae22,Detaf22,Cz22,Dze22,Dzf22,~,~] = getN(SS.A(:,:,2),SS.Bu(:,:,2),SS.Bw(:,:,2),SS.Cy(:,:,2),SS.Dyv(:,:,2),eigL(2),a,b,true,false,wnoise);
[Aol2N,BF2N,Be2N,Bf2N,Ceta2N,Detae2N,Detaf2N,Cz2N,Dze2N,Dzf2N,~,~] = getN(SS.A(:,:,2),SS.Bu(:,:,2),SS.Bw(:,:,2),SS.Cy(:,:,2),SS.Dyv(:,:,2),eigL(N),a,b,true,false,wnoise);


%% LMI Options
LMIOptions = sdpsettings('savesolveroutput',0,'verbose',1,'debug', 1,'solver','sdpt3');
LMIOptions.sdpt3.maxit = 1500;
LMIOptions.sdpt3.steptol= 1e-4;



%% Symbolic decision variables
S       = sdpvar(n.x+2*n.y,n.x+2*n.y,'symmetric');
K_1     = sdpvar(n.u,n.x+n.y,'full');
K1      = [K_1, zeros(n.u,n.y)];
K_2     = sdpvar(n.u,n.x+n.y,'full');
K2      = [K_2, zeros(n.u,n.y)];
G_11    = sdpvar(n.x+n.y,n.x+n.y,'full');
G_12    = sdpvar(n.y,n.y,'full');
G1      = [G_11, zeros(n.x+n.y,n.y); zeros(n.y,n.x+n.y), G_12];
G_21    = sdpvar(n.x+n.y,n.x+n.y,'full');
G_22    = sdpvar(n.y,n.y,'full');
G2      = [G_21, zeros(n.x+n.y,n.y); zeros(n.y,n.x+n.y), G_22];
tt      = sdpvar; % gamma^2
sigmax	= 1/sigma;

%% LMI matrix vt_1 lambda_2 
ACG12 = [Aol12; Ceta12; Cz12];
ACK12 = zeros(size([Aol12; Ceta12; Cz12],1),n.u);
ACK12(1:size(Aol12,1),1:n.u) = BF12;
BD12 = [...
    Be12,    Bf12;
    Detae12,	Detaf12;
    Dze12,   Dzf12];
LMI12 = [blkdiag(G1'+G1-S,eye(size(Be12,2)),tt*eye(size(Bf12,2))),[G1'*ACG12'+K1'*ACK12';BD12'];
    ACG12*G1+ACK12*K1,BD12, blkdiag(S,sigmax*eye(size(Ceta12,1)),eye(size(Cz12,1)))];

%% LMI matrix vt_1 lambda_N
ACG1N = [Aol1N; Ceta1N; Cz1N];
ACK1N = zeros(size([Aol1N; Ceta1N; Cz1N],1),n.u);
ACK1N(1:size(Aol1N,1),1:n.u) = BF1N;
BD1N = [...
    Be1N,     Bf1N;
    Detae1N,  Detaf1N;
    Dze1N,    Dzf1N];
LMI1N = [blkdiag(G1'+G1-S,eye(size(Be1N,2)),tt*eye(size(Bf1N,2))),[G1'*ACG1N'+K1'*ACK1N';BD1N'];
    ACG1N*G1+ACK1N*K1,BD1N, blkdiag(S,sigmax*eye(size(Ceta1N,1)),eye(size(Cz1N,1)))];

%% LMI matrix vt_2 lambda_2
ACG22 = [Aol22; Ceta22; Cz22];
ACK22 = zeros(size([Aol22; Ceta22; Cz22],1),n.u);
ACK22(1:size(Aol22,1),1:n.u) = BF22;
BD22 = [...
    Be22,    Bf22;
    Detae22,	Detaf22;
    Dze22,   Dzf22];
LMI22 = [blkdiag(G2'+G2-S,eye(size(Be22,2)),tt*eye(size(Bf22,2))),[G2'*ACG22'+K2'*ACK22';BD22'];
    ACG22*G2+ACK22*K2,BD22, blkdiag(S,sigmax*eye(size(Ceta22,1)),eye(size(Cz22,1)))];

%% LMI matrix vt_2 lambda_N
ACG2N = [Aol2N; Ceta2N; Cz2N];
ACK2N = zeros(size([Aol2N; Ceta2N; Cz2N],1),n.u);
ACK2N(1:size(Aol2N,1),1:n.u) = BF2N;
BD2N = [...
    Be2N,     Bf2N;
    Detae2N,  Detaf2N;
    Dze2N,    Dzf2N];
LMI2N = [blkdiag(G2'+G2-S,eye(size(Be2N,2)),tt*eye(size(Bf2N,2))),[G2'*ACG2N'+K2'*ACK2N';BD2N'];
    ACG2N*G2+ACK2N*K2,BD2N, blkdiag(S,sigmax*eye(size(Ceta2N,1)),eye(size(Cz2N,1)))];

%% Constraints
LMIConstraints = [tt>=LMImargin,...
    S     >= LMImargin*eye(size(S)),...
    LMI12 >= LMImargin*eye(size(LMI12)),...
    LMI1N >= LMImargin*eye(size(LMI1N)),...
    LMI22 >= LMImargin*eye(size(LMI22)),...
    LMI2N >= LMImargin*eye(size(LMI2N))];

%% Solve LMI
solution  = optimize(LMIConstraints, [], LMIOptions);
fprintf('\n\n***************************************************************************\n\n')
if solution.problem == 0
    fprintf('\nFeasible Solution.\n\n')
else
    solution.info
    yalmiperror(solution.problem)
end

%% Transform matrices

S       = double(S);
G1      = double(G1);
G2      = double(G2);
K1      = double(K1);
K2      = double(K2);
tt      = double(tt);
gamma   = sqrt(tt);
fprintf('Performance threshold: gamma = %d.\n\n',gamma)
sigmax  = double(sigmax);
sigma   = 1/sigmax;
F(:,:,1)        = K1*(inv(G1));
Fx(:,:,1)       = F(:,1:n.x,1);
Fzeta(:,:,1)    = F(:,n.x+1:n.x+n.y,1);
Fz(:,:,1)       = F(:,n.x+n.y+1:end,1);
F(:,:,2)        = K2*(inv(G2));
Fx(:,:,2)       = F(:,1:n.x,2);
Fzeta(:,:,2)    = F(:,n.x+1:n.x+n.y,2);
Fz(:,:,2)       = F(:,n.x+n.y+1:end,2);
[Fx_0,    Fx_1 ]    = split(   Fx(:,:,1),   Fx(:,:,2),vt);
[Fzeta_0, Fzeta_1 ] = split(Fzeta(:,:,1),Fzeta(:,:,2),vt);
[Fz_0,    Fz_1 ]    = split(   Fz(:,:,1),   Fz(:,:,2),vt);

%% Eigenvalue check
% LMIs:
LMI12 = [blkdiag(G1'+G1-S,eye(size(Be12,2)),tt*eye(size(Bf12,2))),[G1'*ACG12'+K1'*ACK12';BD12'];
    ACG12*G1+ACK12*K1,BD12, blkdiag(S,sigmax*eye(size(Ceta12,1)),eye(size(Cz12,1)))];
LMI1N = [blkdiag(G1'+G1-S,eye(size(Be1N,2)),tt*eye(size(Bf1N,2))),[G1'*ACG1N'+K1'*ACK1N';BD1N'];
    ACG1N*G1+ACK1N*K1,BD1N, blkdiag(S,sigmax*eye(size(Ceta1N,1)),eye(size(Cz1N,1)))];
LMI22 = [blkdiag(G2'+G2-S,eye(size(Be22,2)),tt*eye(size(Bf22,2))),[G2'*ACG22'+K2'*ACK22';BD22'];
    ACG22*G2+ACK22*K2,BD22, blkdiag(S,sigmax*eye(size(Ceta22,1)),eye(size(Cz22,1)))];
LMI2N = [blkdiag(G2'+G2-S,eye(size(Be2N,2)),tt*eye(size(Bf2N,2))),[G2'*ACG2N'+K2'*ACK2N';BD2N'];
    ACG2N*G2+ACK2N*K2,BD2N, blkdiag(S,sigmax*eye(size(Ceta2N,1)),eye(size(Cz2N,1)))];

fprintf('***************************************************************************\n\n')
if min([eig(LMI12);eig(LMI1N);eig(LMI22);eig(LMI2N)])>0
    fprintf('LMI was solved successfully!\nSmallest Eigenvalue: %d > 0\nPositive definite Matrix.\n\n',min([eig(LMI12);eig(LMI1N)]))
else
    fprintf('LMI was solved unsuccessfully!\nSmallest Eigenvalue: %d <= 0\nSolution is not positive definite.\n\n',min([eig(LMI12);eig(LMI1N)]))
end
fprintf('***************************************************************************\n\n')

% Closed loop matrices:
Acl2 = Aol12 + BF12*F(:,:,1);
AclN = Aol1N + BF1N*F(:,:,1);
Ac22 = Aol22 + BF22*F(:,:,2);
Ac2N = Aol2N + BF2N*F(:,:,2);
if max(abs([eig(Acl2);eig(AclN);eig(Ac22);eig(Ac2N)])) >= 1
    fprintf('WARNING: The closed loop system is unstable!\nLargest Eigenvalue: |lambda| = %d >= 1.\n',max(abs([eig(Acl2);eig(AclN);eig(Ac22);eig(Ac2N)])))
else
    fprintf('The closed loop system is stable!\nLargest Eigenvalue: |lambda| = %d < 1.\n\n',max(abs([eig(Acl2);eig(AclN);eig(Ac22);eig(Ac2N)])))
end
fprintf('***************************************************************************\n\n')

%% Determine kfinal and t
kfinal = Tfinal/Ts + 1;
t = linspace(0, Tfinal, kfinal);

%% External signal
% Reference
r = rk*ones(size(t));

% Input noise
w = zeros(N*n.w,kfinal);
if wnoise ~= 0
    w = w + wgn(N*n.w,kfinal,10*log10(wnoise)); % in W
    %     %% Uncomment if Noise should be checked
    %     wnoise
    %     mean(var(w))
end
if dist ~= 0
    w(1,find(t == 4):find(t == 4.5)) = w(1,find(t == 4):find(t == 4.5)) + 10*ones(1,find(t == 4.5)-find(t == 4)+1);
end
f = [r;w];

%% Initial states
x           = zeros(N*n.x,kfinal+1);    % states
x(:,1)      = [x01;x02;x03];            % initial states
zeta        = zeros(N*n.y,kfinal+1);    % estimate of the integrated formation error
ztilde      = zeros(N*n.y,kfinal+1);    % observerstates
psi         = [x;zeta;ztilde];          % gathered states
u           = zeros(N*n.u,kfinal);      % input


xhat        = zeros(N*n.x,kfinal+1);    % estimated state
xthat       = zeros(N*n.x,kfinal+1);    % test estimate state
yhat        = zeros(N*n.y,kfinal+1);    % estimated output
ythat       = zeros(N*n.y,kfinal+1);    % test estimate output

et          = zeros(N*n.y,kfinal);      % test estimate error
e           = zeros(N*n.y,kfinal);      % test estimate error
y           = zeros(N*n.y,kfinal);      % intial output
eta         = zeros(N*n.y,kfinal);      % formation error
etahat      = zeros(N*n.y,kfinal);      % estimate of the formation error
etathat     = zeros(N*n.y,kfinal);      % test estimate of the formation error
z           = zeros(N*n.y,kfinal);      % performance output
trigcount   = zeros(N,kfinal);          % amount of triggerings

pos         = zeros(2*N,1);
poshat      = zeros(2*N,1);
posthat     = zeros(2*N,1);

the(:,1)    = zeros(N,1);
thehat      = zeros(N,1);
thethat     = zeros(N,1);
R_the       = eye(2*N);

R_thethat	= zeros(n.y,n.y,N,kfinal);
zetahat     = zeros(N*n.y,kfinal+1);    % estimate of the integrated formation error
zetathat    = zeros(N*n.y,kfinal+1);    % test estimate of the integrated formation error

%% Kronecker Matrices
Q           = eye(N);
Lkron       = kron(L,eye(n.y));
SSkron.Bu	= kron(eye(N),SS.Bu(:,:,nrho));
SSkron.Bw	= kron(eye(N),SS.Bw(:,:,nrho));
SSkron.Cy	= kron(eye(N),SS.Cy(:,:,nrho));
SSkron.Dyv	= kron(eye(N),SS.Dyv(:,:,nrho));

%% k = 1 +
for k = 1:kfinal
    %% Trigger Condition Evaluation
    y(:,k) = SSkron.Cy*x(:,k);
    pos(:,k) = R_the(:,:,k)'*y(:,k);
    et(:,k) = posthat(:,k)-pos(:,k);
    % Trigger condition
    etathat(:,k) = R_the(:,:,k)*Lkron*(r(:,k)-posthat(:,k));
    for ii = 1:N
        if et(n.y*(ii-1)+1:n.y*ii,k)'*et(n.y*(ii-1)+1:n.y*ii,k) <= ETCmargin + sigma*etathat(n.y*(ii-1)+1:n.y*ii,k)'*etathat(n.y*(ii-1)+1:n.y*ii,k)
            switch esti % 0 = ZOH, 1 = OLE, 2 = CLE
                case 0 % ZOH
                    poshat( n.y*(ii-1)+1:n.y*ii,k) = posthat(n.y*(ii-1)+1:n.y*ii,k);
                case 1 % OLE
                    xhat(   n.x*(ii-1)+1:n.x*ii,k) = xthat(  n.x*(ii-1)+1:n.x*ii,k);
                    yhat(   n.y*(ii-1)+1:n.y*ii,k) = ythat(  n.y*(ii-1)+1:n.y*ii,k);
                    thehat( ii,                 k) = thethat(    ii,             k);
                    poshat( n.y*(ii-1)+1:n.y*ii,k) = posthat(n.y*(ii-1)+1:n.y*ii,k);
                case 2 % CLE
                    xhat(   n.x*(ii-1)+1:n.x*ii,k) = xthat(  n.x*(ii-1)+1:n.x*ii,k);
                    yhat(   n.y*(ii-1)+1:n.y*ii,k) = ythat(  n.y*(ii-1)+1:n.y*ii,k);
                    zetahat(n.y*(ii-1)+1:n.y*ii,k) = zetathat(  n.y*(ii-1)+1:n.y*ii,k);
                    thehat( ii,                 k) = thethat(    ii,             k);
                    poshat( n.y*(ii-1)+1:n.y*ii,k) = posthat(n.y*(ii-1)+1:n.y*ii,k);
            end
        else
            trigcount(ii,k) = 1;
            switch esti % 0 = ZOH, 1 = OLE, 2 = CLE
                case 0 % ZOH
                    poshat(n.y*(ii-1)+1:n.y*ii,k) = pos(n.y*(ii-1)+1:n.y*ii,k);
                case 1 % OLE
                    yhat(  n.y*(ii-1)+1:n.y*ii,k) = y(  n.y*(ii-1)+1:n.y*ii,k);
                    xhat(  n.x*(ii-1)+1:n.x*ii,k) = x(  n.x*(ii-1)+1:n.x*ii,k);
                    thehat(ii,                 k) = the(ii,                 k);
                    poshat(n.y*(ii-1)+1:n.y*ii,k) = pos(n.y*(ii-1)+1:n.y*ii,k);
                case 2 % CLE
                    yhat(  n.y*(ii-1)+1:n.y*ii,k) = y(  n.y*(ii-1)+1:n.y*ii,k);
                    xhat(  n.x*(ii-1)+1:n.x*ii,k) = x(  n.x*(ii-1)+1:n.x*ii,k);
                    zetahat(n.y*(ii-1)+1:n.y*ii,k) = zeta(  n.y*(ii-1)+1:n.y*ii,k);
                    thehat(ii,                 k) = the(ii,                 k);
                    poshat(n.y*(ii-1)+1:n.y*ii,k) = pos(n.y*(ii-1)+1:n.y*ii,k);
            end
        end
    end
    
    %% A and F Matrices
    A_help      = zeros(N*n.x);
    Fx_help     = zeros(N*n.u,N*n.x);
    Fzeta_help  = zeros(N*n.u,N*n.y);
    Fz_help     = zeros(N*n.u,N*n.z);
    for ii = 1:N
        A_help(    n.x*(ii-1)+1:n.x*ii,n.x*(ii-1)+1:n.x*ii) = SS.A_0  + SS.A_1  * x(n.x*(ii-1)+4,k);
        Fx_help(   n.u*(ii-1)+1:n.u*ii,n.x*(ii-1)+1:n.x*ii) = Fx_0    + Fx_1    * x(n.x*(ii-1)+4,k);
        Fzeta_help(n.u*(ii-1)+1:n.u*ii,n.y*(ii-1)+1:n.y*ii) = Fzeta_0 + Fzeta_1 * x(n.x*(ii-1)+4,k);
        Fz_help(   n.u*(ii-1)+1:n.u*ii,n.z*(ii-1)+1:n.z*ii) = Fz_0    + Fz_1    * x(n.x*(ii-1)+4,k);
    end
    SSkron.A(:,:,k)    = A_help;
    
    %% Network Dynamics
    eta(:,k) = R_the(:,:,k)*Lkron*(r(:,k)-pos(:,k));
    etahat(:,k) = R_the(:,:,k)*Lkron*(r(:,k)-poshat(:,k));
    u(:,k) = Fx_help*x(:,k) + Fzeta_help*zeta(:,k);% + Fz_help*ztilde(:,k);
    
    % Increment
    x(:,k+1) = SSkron.A(:,:,k)*x(:,k) + SSkron.Bu*u(:,k)  + SSkron.Bw*w(:,k);
    
    zeta(:,k+1) = zeta(:,k)+etahat(:,k);
    ztilde(:,k+1) = - a * ztilde(:,k) + b*(R_the(:,:,k)*r(:,k)-y(:,k));
    
    %% Rotation Matrices
    for ii = 1:N
        the(ii,k+1)=the(ii,k)+(Ts/d)*x(n.x*(ii-1)+4,k);
        R_the(n.y*(ii-1)+1:n.y*ii,n.y*(ii-1)+1:n.y*ii,k+1) = [
            cos(the(ii,k+1))  sin(the(ii,k+1));
            -sin(the(ii,k+1)) cos(the(ii,k+1))];
    end
    
    %% Next Estimate
    switch esti % 0 = ZOH, 1 = OLE, 2 = CLE
        case 0 % ZOH            
            posthat(:,k+1) = poshat(:,k);
        case 1 % OLE
            for ii = 1:N
                xthat(  n.x*(ii-1)+1:n.x*ii,k+1) = (SS.A_0  + (SS.A_1)  * xhat(n.x*(ii-1)+4,k)) *xhat( n.x*(ii-1)+1:n.x*ii,k);
                
                %                 xthat(  n.x*(ii-1)+1:n.x*ii,k+1) = (SS.A_0+SS.Bu(:,:,1)*Fx_0  + (SS.Bu(:,:,1)*Fx_1+SS.A_1)  * xhat(n.x*(ii-1)+4,k)) *xhat( n.x*(ii-1)+1:n.x*ii,k);
                
                ythat(  n.y*(ii-1)+1:n.y*ii,k+1) = SS.Cy(:,:,2)*xthat(n.x*(ii-1)+1:n.x*ii,k+1);
                
                thethat(ii,                 k+1) = thehat(ii,k)+(Ts/d)*xhat(n.x*(ii-1)+4,k);
                
                posthat(n.y*(ii-1)+1:n.y*ii,k+1) = [
                    cos(thethat(ii,k+1))  sin(thethat(ii,k+1));
                    -sin(thethat(ii,k+1)) cos(thethat(ii,k+1))]'*ythat(n.y*(ii-1)+1:n.y*ii,k+1);
            end
        case 2 % CLE
            
            for ii = 1:N
                etahathat(n.y*(ii-1)+1:n.y*ii,k  ) = R_thethat(:,:,ii,k)*kron(Q(ii,:),eye(n.y))*Lkron*(r(:,k)-poshat(:,k));
                Ahat(     :,:,             ii,k  ) = SS.A_0  + SS.A_1  * xhat(n.x*(ii-1)+4,k);
                Fxhat(    :,:,             ii,k  ) = Fx_0    + Fx_1    * xhat(n.x*(ii-1)+4,k);
                Fzetahat( :,:,             ii,k  ) = Fzeta_0 + Fzeta_1 * xhat(n.x*(ii-1)+4,k);
                uhat(     n.u*(ii-1)+1:n.u*ii,k  ) = Fxhat(:,:,ii,k) * xhat( n.x*(ii-1)+1:n.x*ii,k) + Fzetahat(:,:,ii,k)*zetahat(n.y*(ii-1)+1:n.y*ii,k);
                
                zetathat( n.y*(ii-1)+1:n.y*ii,k+1) = zetahat(n.y*(ii-1)+1:n.y*ii,k)+etahathat(n.y*(ii-1)+1:n.y*ii,k);
                xthat(    n.x*(ii-1)+1:n.x*ii,k+1) = Ahat(:,:,ii,k)*xhat( n.x*(ii-1)+1:n.x*ii,k) + SS.Bu(:,:,1)*uhat(n.u*(ii-1)+1:n.u*ii,k);
                ythat(    n.y*(ii-1)+1:n.y*ii,k+1) = SS.Cy(:,:,2)*xthat(n.x*(ii-1)+1:n.x*ii,k+1);
                thethat(                   ii,k+1) = thehat(ii,k)+(Ts/d)*xhat(n.x*(ii-1)+4,k);
                R_thethat(:,:,             ii,k+1) = [
                    cos( thethat(ii,k+1)) sin(thethat(ii,k+1));
                    -sin(thethat(ii,k+1)) cos(thethat(ii,k+1))];
                posthat(  n.y*(ii-1)+1:n.y*ii,k+1) = R_thethat(:,:,ii,k+1)'*ythat(n.y*(ii-1)+1:n.y*ii,k+1);
            end
    end
end
posesti(:,:,esti+1) = pos;
d       = .6;      % 
m       = 1;        % Mass
Iz      = 1;        % 
vt      = 1*[-1,1];   % Scheduling Parameter
SS.A    = zeros(4,4,length(vt));
SS.Bu   = zeros(4,2,length(vt));
SS.Bw   = zeros(4,4,length(vt));
SS.Cy   = zeros(2,4,length(vt));
SS.Dyv  = zeros(2,2,length(vt));

for nrho=1:length(vt)
    SS.A(:,:,nrho) = eye(4) + Ts*[...
        0 vt(nrho)/d 1 0
        -vt(nrho)/d 0 0 1
        0 0 0 0
        0 0 0 0];
    SS.Bu(:,:,nrho) = Ts*[...
        0 0
        0 0
        1/m 0
        0 d/Iz];
    SS.Bw(:,:,nrho) = Ts*eye(4);
    SS.Cy(:,:,nrho)= ... 
        [1 0 0 0
        0 1 0 0];
    SS.Dyv(:,:,nrho) = 10^-3*eye(2);
end


[SS.A_0, SS.A_1 ]= split(SS.A(:,:,1),SS.A(:,:,2),vt);

n.u = 2;
n.v = 2;
n.w = 4;
n.x = 4;
n.y = 2;
n.z = 2;


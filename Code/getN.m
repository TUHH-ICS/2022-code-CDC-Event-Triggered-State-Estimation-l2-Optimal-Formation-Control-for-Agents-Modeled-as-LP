function [Aol,BF,Be,Bf,Ceta,Detae,Detaf,Cz,Dze,Dzf,Axe,Bxe] = getN(A,Bu,Bw,Cy,Dyv,L,a,b,available,v,w)
%sizes
nu = size(Bu,2);
nv = size(Dyv,2);
nw = size(Bw,2);
nx = size(A,1);
ny = size(Cy,1);
%matrix parameters
Az = a*eye(ny);
Bz = b*eye(ny);
if size(L,1) == 1
    L = L*eye(ny);
end
BF1 = Bu;
BF2 = Bu;
BF3 = zeros(2*ny,nu);

A11 = A;
A12 = zeros(nx,nx);
A13 = zeros(nx,2*ny);
A21 = zeros(nx,nx);
A22 = A;
A23 = zeros(nx,2*ny);
A31 = [-L*Cy;-Bz*Cy];
A32 = zeros(2*ny,nx);
A33 = blkdiag(eye(ny),-Az);
Ceta1 = -L*Cy;
Ceta2 = zeros(ny,nx);
Ceta3 = zeros(ny,2*ny);
Cz1 = zeros(ny,nx);
Cz2 = zeros(ny,nx);
Cz3 = [zeros(ny), L];

Be1 = zeros(nx,ny);
Be2 = zeros(nx,ny);
Be3 = [-L;zeros(ny)];
Detae = - L;
Dze = zeros(ny);

Axe = [-Cy,Cy,zeros(ny,2*ny)];
Bxe = [zeros(ny,2*ny)];

Bf1r = zeros(nx,ny);
Bf2r = zeros(nx,ny);
Bf3r = [L;Bz];
Detafr = L;
Dzfr = zeros(ny);
if v == 0
    Bf1v = [];
    Bf2v = [];
    Bf3v = [];
    Detafv = [];
    Dzfv = [];
else
    Bxe = [Bxe, -Dyv];
    
    Bf1v = zeros(nx,nv);
    Bf2v = zeros(nx,nv);
    Bf3v = [-L*Dyv;-Bz*Dyv];
    Detafv = -L*Dyv;
    Dzfv = zeros(ny,nv);
end
if w == 0
    Bf1w = [];
    Bf2w = [];
    Bf3w = [];
    Detafw = [];
    Dzfw = [];
else
    Bxe = [Bxe, zeros(ny,nw)];
    
    Bf1w = Bw;
    Bf2w = zeros(nx,nw);
    Bf3w = zeros(2*ny,nw);
    Detafw = zeros(ny,nw);
    Dzfw = zeros(ny,nw);
end
Bf1 = [Bf1r,Bf1v,Bf1w];
Bf2 = [Bf2r,Bf2v,Bf2w];
Bf3 = [Bf3r,Bf3v,Bf3w];
Detaf = [Detafr,Detafv,Detafw];
Dzf = [Dzfr,Dzfv,Dzfw];


if available
    BF = [BF1;BF3];
    Aol = [A11,A13;A31,A33];
    Be = [Be1;Be3];
    Bf = [Bf1;Bf3];
    Ceta = [Ceta1,Ceta3];
    Cz = [Cz1,Cz3];
else
    BF = [BF1;BF2;BF3];
    Aol = [A11,A12,A13;
        A21,A22,A23;
        A31,A32,A33];
    Be = [Be1;Be2;Be3];
    Ceta = [Ceta1,Ceta2,Ceta3];
    Cz = [Cz1,Cz2,Cz3];
    Bf = [Bf1;Bf2;Bf3];
end


clc
clear
close all

%% physical properties
mu = 1e-2;
rho =1 ;

%% grid
dx = 0.02;
dr = 0.02;
nTubeAxis = 50;
nTubeRadius = 50;
n_points_x = nTubeAxis +1;
n_points_y = nTubeRadius +1;

%% boundary condition
uLeft = 0;
uRight = 0;
uTop = 1;


vBottom = 0;
vTop = 0;

%% initial u & v
u = 0*ones((nTubeAxis+1),nTubeRadius+2); % 11*7 unknown u points
u(:,nTubeRadius+2) = 2;
du = 0*u;

v = 0*ones(nTubeAxis+2,(nTubeRadius+1)); % 40 unknown v points
uGuess = u;
uNew = u;
vGuess = v;
vNew = v;
dv = v;
p = 0*ones(nTubeAxis+2,nTubeRadius+2); % 50 unknown p nodes
pCor = p;
b = pCor;
% cavity_initPressure


%% assign edge point (boundary number)
% cavity_uBoundaryNum;
% cavity_vBoundaryNum;
% cavity_pBoundaryNum;

%% initialize the coeff Mat and source Vec
uMat = zeros(numel(u),numel(u));
uVec = zeros(numel(u),1);
vMat = zeros(numel(v),numel(v));
vVec = zeros(numel(v),1);
pCorMat = zeros(numel(p),numel(p));
pCorVec = zeros(numel(p),1);

pNew = p;
iter = 0;
for exIter = 1:10000
    %% construct uMat and uVec
    for i = 2:n_points_x-1
        for k = 2:n_points_y
            j = 2+n_points_y-k;
            Fe = dr * rho* (u(i,j) + u(i+1,j))/2;
            Fw = dr * rho* (u(i,j) + u(i-1,j))/2;
            Fs = dx * rho* (v(i,j-1) + v(i+1,j-1))/2;
            Fn = dx * rho* (v(i,j) + v(i+1,j))/2;% n or s is a question
            
            De = dr * mu/dx;
            Dw = De;
            Dn = dx * mu/dr;
            Ds = Dn;
            
            du(i,j) = dr/((Dw - Fw/2) + (De+Fe/2) + (Dn  +Fn/2)+ (Ds - Fs/2));
            
            uGuess(i,j) = (Dw + Fw/2)*u(i-1,j) + (De-Fe/2)*u(i+1,j) + (Dn  -Fn/2)*u(i,j+1) + (Ds + Fs/2)*u(i,j-1);
            uGuess(i,j) = uGuess(i,j)/((Dw - Fw/2) + (De+Fe/2) + (Dn  +Fn/2)+ (Ds - Fs/2));
            uGuess(i,j) = uGuess(i,j) + dr*(p(i,j)-p(i+1,j))/((Dw - Fw/2) + (De+Fe/2) + (Dn  +Fn/2)+ (Ds - Fs/2));
        end
    end
    uGuess(:,nTubeRadius+2) = 2-uGuess(:,nTubeRadius+1);
    uGuess(:,1) = -uGuess(:,2);
    uGuess(1,2:nTubeRadius+1) = 0;
    uGuess(nTubeAxis+1,2:nTubeRadius+1) = 0;
    %%
    for i = 2:n_points_x
        for k = 2:n_points_y-1
            j = 2+n_points_y-1-k;
            Fn = rho*dx*(v(i,j) + v(i,j+1))/2;
            Fs = rho*dx*(v(i,j) + v(i,j-1))/2;
            Fe = rho*dr*(u(i,j) + u(i,j+1))/2;
            Fw = rho*dr*(u(i-1,j) + u(i-1,j+1))/2 ;
            Dn = mu*dx/dr;
            Ds = Dn;
            De = mu*dr/dx;
            Dw = De;
            
            vGuess(i,j) = (Dw + Fw/2)*v(i-1,j) + (De -Fe/2)*v(i+1,j) + (Ds + Fs/2)*v(i,j-1) + (Dn - Fn/2)*v(i,j+1);
            vGuess(i,j) = vGuess(i,j)/((Dw - Fw/2) + (De+Fe/2) + (Dn +Fn/2)+ (Ds - Fs/2));
            vGuess(i,j) = vGuess(i,j) + dx*(p(i,j)-p(i,j+1))/((Dw - Fw/2) + (De+Fe/2) + (Dn +Fn/2)+ (Ds - Fs/2));
            
            dv(i,j) = dx/((Dw - Fw/2) + (De+Fe/2) + (Dn +Fn/2)+ (Ds - Fs/2));
        end
    end
    vGuess(1,:) = -vGuess(2,:);
    vGuess(n_points_x+1,:) = -vGuess(n_points_x,:);
    vGuess(2:n_points_x,1) = 0;
    vGuess(2:n_points_x,n_points_y) = 0;

 pCor = 0*pCor;%!!!!!!!!!!!!!!!!!!!!!!!cao 就是这玩意卡了我好久 气死！
    for   i=2:n_points_x
        for  k=2:n_points_y
            j = n_points_y+2-k;
            
            aE = rho*du(i,j)*dr;
            aW = rho*du(i-1,j)*dr;
            
            
            aN = rho*dv(i,j)*dx;
            aS = rho*dv(i,j-1)*dx;
            aP = aE+aW+aS+aN;
            b(i,j) = rho*(dr*uGuess(i-1,j) - dr*uGuess(i,j) + dx*vGuess(i,j-1) - dx*vGuess(i,j) );
            
            
            pCor(i,j) = (aE*pCor(i+1,j) + aW*pCor(i-1,j) + aN*pCor(i,j+1) + aS*pCor(i,j-1) +b(i,j))/aP;

           
            
        end
    end

    
    for i =2:n_points_x
        for j = 2:n_points_y
            
            pNew(i,j) = pNew(i,j) + 0.8*pCor(i,j);
        end
    end
   
    
    pNew(:,n_points_y+1) = pNew(:,n_points_y);
    pNew(:,1) = pNew(:,2);
    pNew(1,:) = pNew(2,:);
    pNew(n_points_x+1,:) = pNew(n_points_x,:);
    
    for i = 2:n_points_x-1
        for j = 2:n_points_y+1
            uNew(i,j) = uGuess(i,j) + 0.8*du(i,j)*(pCor(i,j)-pCor(i+1,j));
            
        end
    end
    uNew(:,nTubeRadius+2) = 2-uNew(:,nTubeRadius+1);
    uNew(:,1) = -uNew(:,2);
    uNew(1,2:nTubeRadius+1) = 0;
    uNew(nTubeAxis+1,2:nTubeRadius+1) = 0;
    
    for i = 2:n_points_x
        for j = 2:n_points_y-1
            vNew(i,j) = vGuess(i,j) + 0.8*dv(i,j)*(pCor(i,j)-pCor(i,j+1));
        end
    end
    vNew(1,:) = -vNew(2,:);% left
    vNew(nTubeAxis+2,:) = -vNew(nTubeAxis+1,:);%right
    vNew(2:nTubeAxis+1,1) = 0;
    vNew(2:nTubeAxis+1,nTubeRadius+1) = 0;
    
    
    error = 0;
    for i = 2:n_points_x
        for j = 2:n_points_y
            error = error + abs(b(i,j));
        end
    end
    
    u = uNew;
    v = vNew;
    p = pNew;
    iter = iter +1;
    if mod(iter,100)==0
        contour(u(2:n_points_x,2:n_points_y));
        drawnow
        colorbar
        error
        %     pause
    end
    if error < 1e-5
        break
    end
    
end

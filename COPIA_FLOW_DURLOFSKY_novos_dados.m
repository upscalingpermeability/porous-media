%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local-global upscaling based on Y. Chen, L.J. Durlofsky, M. Gerritsen, 
%% and X.H. Wen. A coupled local–global upscaling approach for simulating 
%% flow in highly heterogeneous formations. Advances in Water Resources, 
%% 26(10):1041–1060, 2003.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author...: Marcio Borges
%% Date.....: 22/10/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
%addpath ~/Dropbox/mrst-2021b/
addpath ./mrst-2021b/
startup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
   require upscaling coarsegrid incomp mimetic mpfa
catch %#ok<CTCH>
   mrstModule add upscaling coarsegrid incomp mimetic mpfa
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart  = tic;
verbose = false;
color   = 'none';
color2  = 'black';
%color2  = 'none';
lim     = [0 0];
printa  = 1;
prt     = 1; % print for simulation
met     = 'tpfa'; %Discretization method used in computation.
% name = 'gauss_10x10_local_';
%name = 'gauss_6x22_local';
% name = 'campos_10x10_local_teste';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lx  = 120.0;
% Ly  = 440.0;
% Lz  = 1.0;
Lx  = 50.0;
Ly  = 50.0;
Lz  = 1.0;

%% fine grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nx  = 60;
% ny  = 220;
% nz  = 1;
nx  = 50;
ny  = 50;
nz  = 1;
%% coarse grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nnx = 6;
% nny = 22;
% nnz = 1;
nnx = 10;
nny = 10;
nnz = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home   = './campos_artigo';
homef  = './figuras_artigo/';
nini   = 0;
nfim   = 4999;
% filenx = 'fields/analitico_512x512x20_256x256x1_x_';
% fileny = 'fields/analitico_512x512x20_256x256x1_y_';
% filenz = 'fields/analitico_512x512x20_256x256x1_z_';

% filenx = 'campos_artigo5000/campos_';
% fileny = 'campos_artigo5000/campos_';
% filenz = 'campos_artigo5000/campos_';


% filenx = 'fields/e512x512x20_32x32x4_l32x32x8_';
% fileny = 'fields/e512x512x20_32x32x4_l32x32x8_';
% filenz = 'fields/e512x512x20_32x32x4_l32x32x8_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namein    = [home '/dataLG_' num2str(nx, '%d') 'x' num2str(ny, '%d') 'x'...
    num2str(nz, '%d') '_' num2str(nfim-nini+1,'%d') '.mat'];
nameout   = [home '/dataLG_' num2str(nnx,'%d') 'x' num2str(nny,'%d') 'x'...
    num2str(nnz,'%d') '_' num2str(nfim-nini+1,'%d') '.mat'];
fileIDin  = fopen(namein,'w');
fileIDout = fopen(nameout,'w');
fprintf('\n ===============================================================')
fprintf('\n Output Files:\n%s\n%s',namein,nameout);
fprintf('\n ===============================================================\n')
%% fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho  = 1.0
beta = 1.0*milli() * darcy();
raio = 0.5;
%NMAX = 10; % maximum number of iteration
NMAX = 20; % maximum number of iteration
TOL  = 1.0e-4;
%name = 'local_global_full';
et   = 0.5;
fulltensor = false; % if == true => use the MPFA version and retorn a full perm. tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dim, nD, fine_grid, coarse_grid, dims] = preproc(Lx,Ly,Lz,nx,ny,nz,...
    nnx,nny,nnz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fine Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G     = cartGrid(fine_grid, dims);
G     = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coarse grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_ups = cartGrid(coarse_grid, dims);
G_ups = computeGeometry(G_ups);
p     = partitionUI(G, coarse_grid);
p     = processPartition(G, p);
CG    = generateCoarseGrid(G, p, 'Verbose', verbose);
CG    = coarsenGeometry(CG);
%% Main loop ==============================================================
id = [1; 2; 3];
if fulltensor, id = [1; 4; 6]; end
m   = 0;
fat = 1.0;
beta= beta * fat;
for n = nini:nfim
name = ['campos_10x10_local_teste_' num2str(n+1)];
%  
%     load('saida.mat');
%     K=input(:,n+1);
%      load('data.mat');
%     K=MATRIZ(:,n+1);

%     load('dados50_artigo.mat')
%     K=z(:,n+1);
    load('dados5000_artigo.mat')
    K=z(:,n+1);
    rock.perm    = beta * exp(rho * K);
    rockUps.perm = upscalePerm(G, CG, rock, 'Verbose', verbose, 'method', met); 
    if dim == 2
        in_fields    = [K(:,1); K(:,1)];
        out_fields   = [reverseKlog(rockUps.perm(:,1),beta,rho);...
            reverseKlog(rockUps.perm(:,2),beta,rho)];
    else
        in_fields    = [K(:,1); K(:,2); K(:,3)];
        out_fields   = [reverseKlog(rockUps.perm(:,1),beta,rho);...
            reverseKlog(rockUps.perm(:,2),beta,rho);...
            reverseKlog(rockUps.perm(:,3),beta,rho)];
    end
    fwrite(fileIDin,in_fields,'single');
    fwrite(fileIDout,out_fields,'single');
    save_field1(Lx,Ly,Lz,nnx,nny,nnz,1,...
        rockUps.perm,rho,beta,id,n,home,name,dim,nD,prt,'Y');
    m = m+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
telapsed = toc(tstart);
Totaltime= seconds(telapsed);
Totaltime.Format = 'hh:mm:ss';
fprintf('\n =======================================================\n')
fprintf(' Total time elapsed.......: %s\n',char(Totaltime))
fprintf(' =======================================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printa <= 2
    id = [1;2;3];
    if fulltensor, id = [1; 4; 6]; end
    if dim == 2
        %vw  = [90 90];
        vw  = [0 90];
    else
        vw  = [-35 20];
%         vw  = [0 90];
    end
    
%     rockUps.perm = [rockUps.perm(:,id(1)) rockUps.perm(:,id(2))...
%         rockUps.perm(:,id(3))];

    lim = [min(min(reverseKlog(rock.perm,beta(1:1),rho(1:1))))...
        max(max(reverseKlog(rock.perm,beta(1:1),rho(1:1))))];
    
%     lim = [min(min(reverseKlog(rockUps.perm(:,2),beta(1:dim),rho(1:dim))))...
%         max(max(reverseKlog(rockUps.perm(:,2),beta(1:dim),rho(1:dim))))]

    G_ups = cartGrid(coarse_grid, dims);
    G_ups = computeGeometry(G_ups);
    
    if dim == 2
        figure_upscaling22D(rock,rockUps,G,G_ups,beta,rho,color,...
            color2,lim,vw);
    else
        figure_upscaling(rock,rockUps,G,G_ups,beta,rho,color,...
            color2,lim,vw);
    end
    base=[homef 'upscaling_' name];
    set(gcf,'PaperPositionMode','auto');
    print('-dpng','-r600', base);
end
% rockUps.perm
mean(rockUps.perm)
%clear





% perm{1} = reverseKlog(rock.perm(:,1),beta,rho);
% perm{2} = reverseKlog(rockUps.perm(:,1),beta,rho);
% perm{3} = reverseKlog(rockUps.perm(:,2),beta,rho);

% 
% for i=1:3
%    subplot(1,3,i);
%    [h,ax] = colorbarHist(perm{i},[-18 -10],'South');
%    %set(h,'XTick',-16:2:-12,'XTickLabel',{'.1', '10', '1000'});
% end


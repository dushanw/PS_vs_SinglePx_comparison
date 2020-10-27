% 20201025 by Dushan N Wadduwage
% test for PS vs. SinglePixel

clc; clear all; close all
addpath('_datasets');

p           = 1/4;
Nx_list     = [8 16 32];
exp_name    = 'mousebrain';                              % {beads, minist,'mousebrain'} mousebrain downloaded from: https://connectivity.brain-map.org/projection/experiment/siv/515418047?imageId=515418502&imageType=TWO_PHOTON,SEGMENTATION&initImage=TWO_PHOTON&x=20545&y=9347&z=6
exp_type    = 'eq-time';                            % {'eq-time','eq-photons'}
N_reps      = 100;                                  % #noise-loop-repetitions
N_beads     = 5;                                    % used only for beads

res_dir     = ['./_results/' sprintf('/p=%.3f/%s/',p,datetime)];
mkdir(res_dir);

for ii=1:length(Nx_list)

  %% pram init
  Nx          = Nx_list(ii);
  Ny          = Nx;
  compression = 1;                                    % compression >= 1
  M           = Nx*Ny/compression;                    % number of measurements

  alpha       = 2/(M*p);                              % from Eq 6

  %% test objects
  switch exp_name
    case 'minist'
      load ./_datasets/minist.mat                                                         % MINIST data set
      fig_name  = sprintf('%s_%s_NxNy-%dx%d_M-%d_',exp_type,exp_name,Nx,Ny,M);
    case 'beads'
      XTest     = f_genobj_beads3D(Ny*4,Nx*4,50,N_beads);                                 % synth beads
      XTest     = sum(XTest,3);
      fig_name  = sprintf('%s_%s-%d_NxNy-%dx%d_M-%d_',exp_type,exp_name,N_beads,Nx,Ny,M);
    case 'mousebrain'
      XTest     = imread('./_datasets/515418047_40.jpg');                         % mouseBrain downloaded
      XTest     = mean(single(XTest),3);
      XTest     = XTest(350:450 ,600:700);
      fig_name  = sprintf('%s_%s_NxNy-%dx%d_M-%d_',exp_type,exp_name,Nx,Ny,M);
  end
  XTest       = imresize(XTest,[Ny Nx]);
  X           = XTest(:,:,1,1);
  X           = X-min(X(:));
  X           = X/sum(X(:));

  %% measurement matrix
  A           = alpha*(rand(M,Nx*Ny) < p);

  %% image
  switch exp_type
    case 'eq-photons'
      mux_factor  = 1;
      x_label     = 'Normalized #photons [AU]';
    case 'eq-time'
      mux_factor  = Nx*Ny;
      x_label     = 'Normalized time [AU]';
  end

  tic 
  Xhat_I = [];
  for i=1:15
    i
    t_normed(i)     = 10^i;
    totPhotons_ps(i)= 10^i;
    totPhotons_spx  = 10^i * mux_factor;

    Y_spx           = A*X(:) * totPhotons_spx;
    Y_ps            =   X    * totPhotons_ps(i); 
    for j = 1:N_reps
      % SPX
      Yhat_spx      = poissrnd(Y_spx);                      % add noise to spx measurement    
      Xhat_spx      = A\Yhat_spx;                           % reconstruct using mse
      Xhat_spx      = reshape(Xhat_spx,[Ny Nx]);             
      ssim_spx(i,j) = ssim(rescale(Xhat_spx),rescale(X));   % ssim_spx  
      % PS
      Xhat_ps       = poissrnd(Y_ps);                       % point scan image
      ssim_ps (i,j) = ssim(rescale(Xhat_ps),rescale(X));    % ssim_ps  
    end
    Xhat_I          = cat(2,Xhat_I,cat(1,rescale(Xhat_spx),rescale(Xhat_ps)));
    Xhats_spx{i,j}  = Xhat_spx;
    Xhats_ps {i,j}  = Xhat_ps;
  end
  toc
  
  %% figures
  errorbar(totPhotons_ps,mean(ssim_spx,2),std(ssim_spx,0,2),'LineWidth',2);hold on
  errorbar(totPhotons_ps,mean(ssim_ps ,2),std(ssim_ps ,0,2),'LineWidth',2);hold off
  set(gca,'XScale','log')
  xticks(totPhotons_ps(1:2:end))
  ylim([0 1.3])
  ylabel('SSIM value [AU]')
  xlabel(x_label)
  legend({'Single pixel','Point scan'})
  set(gca,'fontsize',24);
  saveas(gcf,sprintf('%s%s_ssim.png',res_dir,fig_name))

  imagesc(Xhat_I);axis image;axis off
  saveas(gcf,sprintf('%s%s_rep-Xhat-spix-ps.png',res_dir,fig_name))

  close all

  %% run results to save    
  run_results_forNx{i}.Xhats_spx      = Xhats_spx;
  run_results_forNx{i}.Xhats_ps       = Xhats_ps;
  run_results_forNx{i}.ssim_spx       = ssim_spx;
  run_results_forNx{i}.ssim_ps        = ssim_ps;
  run_results_forNx{i}.totPhotons_ps  = totPhotons_ps;
  run_results_forNx{i}.Nx             = Nx;
  run_results_forNx{i}.M              = M;
  run_results_forNx{i}.p              = p;

end

save([res_dir exp_name '-' exp_type '_run-results.mat'],'run_results_forNx');

%% Retrieve saved results and plot 

load([res_dir exp_name '-' exp_type '_run-results.mat'],'run_results_forNx');% <replace with the right name>










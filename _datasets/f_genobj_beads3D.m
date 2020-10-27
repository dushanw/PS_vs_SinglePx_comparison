
function X0 = f_genobj_beads3D(Ny,Nx,Nz,N_beads)
  
  r_mean    = round(min([Ny,Nx,Nz])/10);
  r_range   = r_mean-round(r_mean/5):r_mean+round(r_mean/5);  
  
  if isempty(N_beads)
    N_beads   = round((Ny*Nx*Nz/max(r_range)^3)/100);             % 1/100 of grid maximum packing
  end
    
  centroids = [randi([max(r_range), Ny-max(r_range)],N_beads,1),...
               randi([max(r_range), Nx-max(r_range)],N_beads,1),...
               randi([max(r_range), Nz-max(r_range)],N_beads,1)];
             
  radii     = randi([min(r_range), max(r_range)],N_beads,1);
  norm_int  = normrnd(1,0.1,N_beads,1);                           % normalized intensity
  
  X0        = single(zeros(Ny,Nx,Nz));
  
  [Y,X,Z]   = meshgrid(1:Ny,1:Nx,1:Nz);
    
  for i=1:N_beads
    i;
    beads_fll = (Y - centroids(i,1)).^2 + (X - centroids(i,2)).^2 + (Z - centroids(i,3)).^2 < radii(i)^2;
    X0        = X0 + beads_fll*norm_int(i);
  end

  % volshow(X0)  
end
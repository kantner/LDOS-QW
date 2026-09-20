clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical constants
  par.m0   = 9.1093837139 * 1E-31;
  par.e0   = 1.602176634  * 1E-19;
  par.hbar = 1.054571817  * 1E-34;

% units
  par.nm  = 1E-9;
  par.um  = 1E-6;
  par.eV  = par.e0;
  par.meV = 1E-3 * par.eV;

% grid  
  par.L = 500 * par.nm;
  par.N = 1001;
  par.dz = par.L/(par.N-1);

  par.z = [0:par.N-1]'*par.dz - par.L/2;

% material parameters
  par.m = 0.067 * par.m0;
  
% quantum well
  par.w = 10 * par.nm; % QW width

  par.sigma_u = 0.01 * 0.25 * par.nm; % upper interface width
  par.sigma_l = 0.01 * 0.25 * par.nm;   % lower interface width

  par.U0 = 0.2 * par.eV; % QW depth

% QW confinement potential
  sigmoid = @(z) 1./(exp(-z) + 1);

  par.U_QW =  par.U0 * ( sigmoid((par.z-par.w/2)/par.sigma_u) -sigmoid((par.z+par.w/2)/par.sigma_l) );

% electric field
  par.F = 0.0 * 1E6; %(unit: V/m)

% full confinement potential
  par.U = par.U_QW - par.e0 * par.F * par.z;

% energy range for LDOS computation
  par.E_range = linspace(min(par.U), 0.2*par.eV, 10001);
  par.dE      = mean(diff(par.E_range));

  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute LDOS for bound states
  [LDOS_bound, psi_bound, E_bound] = compute_LDOS_bound(par);

% compute LDOS for propagating states    
  [LDOS_prop_l2r, LDOS_prop_r2l] = compute_LDOS_prop(par);

% compute full LDOS
  LDOS_tot = LDOS_prop_r2l + LDOS_prop_l2r + LDOS_bound;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot bound states
  figure(1); clf; hold on;
  plot(par.z/par.nm, par.U/par.eV, 'k-','LineWidth',2)
  box on
  xlabel('z (nm)')
  ylabel('energy (eV)')
  for iE = 1 : length(E_bound)
    plot(par.z/par.nm,E_bound(iE)/par.eV + abs(psi_bound(:,iE)).^2,'r-')
  end
  title('bound states')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot LDOS
  figure(2);clf; hold on;

    DOS_scale =  1/(par.eV * par.nm^3);

    z_level   = max(max(LDOS_tot/DOS_scale)); % level line for mesoscopic potential


  subplot(2,2,1); hold on;  
    surf(par.z/par.nm, par.E_range/par.eV, LDOS_bound'/DOS_scale)
    plot3(par.z/par.nm, par.U/par.eV,z_level * ones(par.N,1),'w--','LineWidth',2)
    shading interp
    box on
    xlabel('z (nm)')
    ylabel('energy (eV)')
    colorbar
    xlim([min(par.z),max(par.z)]/par.nm)
    ylim([min(par.E_range),max(par.E_range)]/par.eV)
    title('LDOS bound states')

  subplot(2,2,2); hold on;
    surf(par.z/par.nm, par.E_range/par.eV, LDOS_prop_l2r'/DOS_scale)
    plot3(par.z/par.nm, par.U/par.eV,z_level * ones(par.N,1),'w-','LineWidth',2)
    shading interp
    box on
    xlabel('z (nm)')
    ylabel('energy (eV)')
    colorbar
    xlim([min(par.z),max(par.z)]/par.nm)
    ylim([min(par.E_range),max(par.E_range)]/par.eV)
    title('LDOS propagating states (incoming from left)')

  subplot(2,2,3); hold on;
    surf(par.z/par.nm, par.E_range/par.eV, LDOS_prop_r2l'/DOS_scale)
    plot3(par.z/par.nm, par.U/par.eV,z_level * ones(par.N,1),'w-','LineWidth',2)
    shading interp
    box on
    xlabel('z (nm)')
    ylabel('energy (eV)')
    colorbar
    xlim([min(par.z),max(par.z)]/par.nm)
    ylim([min(par.E_range),max(par.E_range)]/par.eV)
    title('LDOS propagating states (incoming from right)')

  subplot(2,2,4); hold on;
    surf(par.z/par.nm, par.E_range/par.eV, LDOS_tot'/DOS_scale)
    plot3(par.z/par.nm, par.U/par.eV,z_level * ones(par.N,1),'w-','LineWidth',2)
    shading interp
    box on
    xlabel('z (nm)')
    ylabel('energy (eV)')
    colorbar
    xlim([min(par.z),max(par.z)]/par.nm)
    ylim([min(par.E_range),max(par.E_range)]/par.eV)
    title('full LDOS')

  sgtitle('local density of states')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot DOS  
  figure(3);clf;hold on;


% bulk DOS for reference
  DOS_bulk = @(E) 1/(2*pi^2) * (2*par.m/(par.hbar*par.hbar))^1.5 * sqrt(max(E,0));

% slice at left facet + reference
%{
  plot(par.E_range/par.eV, LDOS_tot(1,:)/DOS_scale,'r','LineWidth',2,'DisplayName','LDOS (left facet)')
  plot(par.E_range/par.eV, DOS_bulk(par.E_range-par.U(1))  /DOS_scale,'r--','LineWidth', 2,'DisplayName','DOS bulk (left)')  
%}

% slice at mid barrier
  z0  = 0.25*(par.L+par.w);
  idx = min(find(par.z>=z0));

  plot(par.E_range/par.eV, LDOS_tot(idx,:)/DOS_scale,'m','LineWidth',2,'DisplayName','LDOS (barrier)')
  plot(par.E_range/par.eV, DOS_bulk(par.E_range-par.U(idx))  /DOS_scale,'m--','LineWidth', 2,'DisplayName','DOS bulk (barrier)') 


% slice at right facet + reference
%
  plot(par.E_range/par.eV, LDOS_tot(end,:)/DOS_scale,'b','LineWidth',2,'DisplayName','LDOS (right facet)')
  plot(par.E_range/par.eV, DOS_bulk(par.E_range-par.U(end))/DOS_scale,'b--','LineWidth', 2,'DisplayName','DOS bulk (right)')
%}

% average over QW
  indicator = -par.U_QW/par.U0;
  DOS_QW    = trapz(par.dz,LDOS_tot.*indicator,1)/trapz(par.dz,indicator);

  plot(par.E_range/par.eV, DOS_QW/DOS_scale,'g-','LineWidth', 2,'DisplayName','QW averaged')
  plot(par.E_range/par.eV, DOS_bulk(par.E_range-par.U(round(par.N/2)) )/DOS_scale,'g--','LineWidth', 2,'DisplayName','DOS bulk (center)')



  legend()
  xlabel('energy (eV)')
  ylabel('DOS (eV^{-1} nm^{-3})')
  box on
  
 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E,psi] = bound_states(par)
% compute bound states

% internal energy rescaling
  energy_scale = par.eV;

% internal energy shift (set confinement potential minimum to zero)
  minU = min(par.U);

% kinetic energy
  prefactor = (par.hbar/par.dz)^2 *1/(2*par.m * energy_scale) ;
  T = -prefactor* diag(sparse(ones(par.N-1,1)),+1) ... % upper diagonal
      -prefactor* diag(sparse(ones(par.N-1,1)),-1) ...% lower diagonal
      +2*prefactor * diag(sparse(ones(par.N,1))); % diagonal

% potential energy
  U = diag(sparse(par.U - minU))/energy_scale;

% Hamiltonian
  H = T + U;

% solve
  neigs   = 10;
  [psi,E] = eigs(H,neigs,'smallestabs');

% rescale
  E = diag(E) * energy_scale + minU; 

% keep only states below continuum band edge
  E_max = min([par.U(1),par.U(end)]);

  idx   = E<E_max;

  E   = E(idx);
  psi = psi(:,idx);

end


function [psi] = propagating_wave(E,direction,par)
% solve scattering problem for propagating waves in energetic continuum
%
% direction: + 1 = incoming wave from left
%            - 1 = incoming wave from right

% internal energy rescaling
  energy_scale = par.eV;

% kinetic energy
  prefactor = (par.hbar/par.dz)^2 *1/(2*par.m*energy_scale) ;
  T = -prefactor* diag(sparse(ones(par.N-1,1)),+1) ... % upper diagonal
      -prefactor* diag(sparse(ones(par.N-1,1)),-1) ...% lower diagonal
      +2*prefactor * diag(sparse(ones(par.N,1))); % diagonal

% potential energy + injection energy level
  D = diag(sparse((par.U - E)/energy_scale));

% BCs (overwrite first and last line of kinetic energy operator)
  RHS = spalloc(par.N,1,1);

  switch(direction)
    case +1 % incoming from left

      k = 2/par.dz * asin(sqrt(0.5*(E-par.U(1))*par.m*(par.dz/par.hbar)^2));

      T(1,1)         = +prefactor * (2 - exp(1i*k*par.dz));
      T(par.N,par.N) = +prefactor * (2 - exp(1i*k*par.dz));
      RHS(1)         = -prefactor * exp(1i*k*par.z(1)) * 2*1i*sin(k*par.dz);

    case -1 % incoming from right

      k = 2/par.dz * asin(sqrt(0.5*(E-par.U(end))*par.m*(par.dz/par.hbar)^2));
      
      T(1,1)         = +prefactor * (2 - exp(1i*k*par.dz));
      T(par.N,par.N) = +prefactor * (2 - exp(1i*k*par.dz));
      RHS(par.N)     = -prefactor * exp(-1i*k*par.z(end)) * 2*1i*sin(k*par.dz);
    otherwise
      error('invalid option')

  end

% build matrix M = H - E
  M = T + D;

% solve M * psi = RHS
  psi = M\RHS;

% check for errors
  if sum(isnan(psi)) > 0
    error('psi contains NaN values')    
  end
 
end


function [LDOS_bound, psi_bound, E_bound] = compute_LDOS_bound(par)

% compute bound states  
  [E_bound,psi_bound] = bound_states(par);

% allocate memory
  LDOS_bound    = zeros(par.N, length(par.E_range));

% compute LDOS for bound states  
  prefactor = par.m/(pi * par.hbar^2);

% sweep over bound states  
  for n = 1 : length(E_bound)  
      En = E_bound(n);
        
    % add states if E >= E_n
      idx = par.E_range >= En;

    % renormalize (physical scaling)
      psi = psi_bound(:,n)/sqrt(trapz(par.z, abs(psi_bound(:,n)).^2 ));

      LDOS_bound(:,idx) = LDOS_bound(:,idx) ...
                          + prefactor * abs(psi).^2  * ones(1,sum(idx));
    
  end

end


function [LDOS_prop_l2r, LDOS_prop_r2l] = compute_LDOS_prop(par)
% compute LDOS for propagating states

% allocate memory
  LDOS_prop_l2r = zeros(par.N, length(par.E_range));
  LDOS_prop_r2l = zeros(par.N, length(par.E_range));

  prefactor = par.m/(pi * par.hbar^2);
  

  k_discrete = @(E) 2/par.dz * asin(sqrt(0.5*E*par.m*(par.dz/par.hbar)^2));


  for iE = 2 : length(par.E_range)
    E1 = par.E_range(iE-1);
    E2 = par.E_range(iE);

    % left to right
    if E2 > par.U(1)

      % clip the lower interval end at the continuum threshold
      E1c = max(E1, par.U(1));

      %k1 = sqrt(2*par.m*(E1c-par.U(1)))/par.hbar;
      %k2 = sqrt(2*par.m*(E2 -par.U(1)))/par.hbar;

      k1 = k_discrete(E1c-par.U(1));
      k2 = k_discrete(E2 -par.U(1));

      k_mid   = 0.5*(k1+k2);
      %eps_mid = par.U(1) + (par.hbar * k_mid)^2 /(2*par.m);
      eps_mid = par.U(1) + 2/par.m * (par.hbar/par.dz)^2 * sin( k_mid * par.dz/2)^2;

      [psi] = propagating_wave(eps_mid,+1,par);    

      weight = (k2-k1)/(2*pi);

      LDOS_prop_l2r(:,iE:end) = LDOS_prop_l2r(:,iE:end) ...
                                + prefactor * weight * abs(psi).^2;
    end
  
    % right to left
    if E2 > par.U(end)
      % clip the upper interval end at the continuum threshold
      E1c = max(E1, par.U(end));

      %k1 = sqrt(2*par.m*(E1c-par.U(end)))/par.hbar
      %k2 = sqrt(2*par.m*(E2 -par.U(end)))/par.hbar

      k1 = k_discrete(E1c-par.U(end));
      k2 = k_discrete(E2 -par.U(end));

      k_mid   = 0.5*(k1+k2);
      %eps_mid = par.U(end) + (par.hbar * k_mid)^2 /(2*par.m);
      eps_mid = par.U(end) + 2/par.m * (par.hbar/par.dz)^2 * sin( k_mid * par.dz/2)^2;

      [psi] = propagating_wave(eps_mid,-1,par);    

      weight = (k2-k1)/(2*pi);

      LDOS_prop_r2l(:,iE:end) = LDOS_prop_r2l(:,iE:end) ...
                                + prefactor * weight * abs(psi).^2;

    end

  end
end
% To produce the axisymmetric shapes shown in Fig 2C, 2F
% (area 4790 um^2, six rim perimeters; each shape is hung from its rim point)
Perimeter = [220, 180, 140, 100, 60, 20]; 
Areas = 4790*ones(length(Perimeter),1);

A_total    = 1;   % desired area
 
for i = 1:length(Areas)

    hole_rad = Perimeter(i)/(2*pi); 
    SA = Areas(i); 
    L0     = sqrt(SA);

   
    R1_desired = hole_rad/L0;   % right-end radius
    n_mesh     = 100;    % initial mesh points

    sol = findSurfaceBVP2(R1_desired,A_total,n_mesh);
   % E_bend(i) = int_bend_energy(sol);

    % Evaluate the solution on a finer grid if needed:
    newpts = 1000;
    xFine = linspace(0,1,newpts);
    yFine = deval(sol,xFine);   % yFine is 6×length(xFine)
    L = sol.parameters(1);

    % Evaluate solution
    psi = yFine(1,:);           % a(x) = psi(s)
    psip = yFine(2,:)./(L0*L);
    s = L * xFine*L0;          % since ds/dx = L

    r = yFine(3,:)*L0;  % 3rd row is r(x)
    z = yFine(4,:)*L0;  % 4th row is z(x)

    rdata = [z z];
    zdata = [-r r];
    [x_rot, y_rot, phi, xcm_rot, ycm_rot] =  hang_curve([z z], [-r r], 2*newpts); %this must be 2*newpts because this the point we rotate about
    
 
  
   figure; axis equal
   hold on
   box on 
   scatter(x_rot,y_rot-y_rot(end))
   hold off
   % filename = sprintf('Shape_%d.png', i);
   % print(gcf, filename, '-dpng', '-r600'); % gcf = current figure, 300 dpi

    
end

%% Mean curvature H(r) of the case-0 shapes -> Cond1.mat ... Cond6.mat
% Same six shapes as case 0 in findSurfaceBVP_backup.m. Each file stores
% X = [r, H] (column vectors): r in um, only r > 0.5 um as in the original
% Cond files, and the mean curvature H = (dpsi/ds + sin(psi)/r)/2 in 1/um.
% To write the files, uncomment the outdir and save lines below; they go to
% the folder Cond_files so the existing Cond*.mat in this folder are not
% overwritten.
Perimeter = [186, 184, 184, 177, 166, 154];      % rim perimeters (um)
Areas     = [2911, 2925, 2988, 3056, 3075, 3086]; % membrane areas (um^2)
%outdir = 'Cond_files';
%if ~exist(outdir, 'dir'), mkdir(outdir); end

figure; hold on; box on
for i = 1:length(Areas)
    L0  = sqrt(Areas(i));                                % length scale (um)
    sol = findSurfaceBVP2(Perimeter(i)/(2*pi)/L0, 1, 100); % rim radius / L0, area 1, 100 mesh pts

    % BVP solution on 1000 points, converted to um
    xFine = linspace(0,1,1000);
    yFine = deval(sol, xFine);
    L     = sol.parameters(1);
    psi   = yFine(1,:);            % tangent angle psi
    psip  = yFine(2,:)./(L0*L);    % dpsi/ds (1/um): BVP stores dpsi/dx, and ds/dx = L*L0
    s     = L*xFine*L0;            % arc length (um)
    r     = yFine(3,:)*L0;         % radius (um)

    H = axisym_mean_curvature(s, psi, psip);
    X = [r(r>0.5)', H(r>0.5)'];
    %save(fullfile(outdir, sprintf('Cond%d.mat', i)), 'X');

    plot(X(:,1), X(:,2), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('P = %d, A = %d', Perimeter(i), Areas(i)));
end
xlabel('$r\ (\mathrm{\mu m})$', 'Interpreter', 'latex');
ylabel('$H\ (\mathrm{\mu m^{-1}})$', 'Interpreter', 'latex');
legend('Location', 'northeast')




%% Continuation sweep in the rim radius R (dimensionless, area = 1)
% Solve the shape BVP at r_init, then step R up to r_max (close to the flat
% disk, R = 1/sqrt(pi) = 0.5642) and down to r_min (nearly closed), each
% time starting from the previous solution.
A      = 1;
r_init = 0.05;
r_min  = 0.0001;       % smallest rim radius (backward sweep)
r_max  = 0.5640;       % largest rim radius, just below the flat disk 1/sqrt(pi)
r_step = 0.001;
n_mesh = 100;          % initial mesh points for the first solve

% Initial guess/solution at r_init (can be solinit or previous bvp4c solution)
sol_prev = findSurfaceBVP2(r_init, A, n_mesh);

% bvp4c default tolerances (RelTol = 1e-3). Checked 2026-09-30: with
% RelTol = AbsTol = 1e-6 the plotted energies change by < 0.1 kBT and the
% critical rhog by at most one rhogvec step; only shapes with L < 0.5 um
% change noticeably.
opts = bvpset('NMax', 1e5);

sweep_sols = cell(0,2);

r = r_init;
while r < r_max - 1e-12
    target = min(r + r_step, r_max);

    % Try full step first
    bc = @(ya,yb,p) bcfun(ya, yb, p, target, A);
    try
        sol = bvp4c(@odefun, bc, sol_prev, opts);
        success = true;
    catch
        success = false;
    end

    % If bvp4c threw an error, retry with up to 6 halvings of the step
    % (bvp4c warnings are not caught here)
    if ~success
        step = (target - r)/2;
        ok = false;
        for k = 1:6
            candidate = r + step;
            bc_cand = @(ya,yb,p) bcfun(ya, yb, p, candidate, A);
            try
                sol = bvp4c(@odefun, bc_cand, sol_prev, opts);
                ok = true; target = candidate;
                break
            catch
                step = step/2;
            end
        end
        if ~ok
            fprintf('stalled near r=%.4f\n', r);
            break
        end
    end

    % accept the step
    sol_prev = sol;
    r = target;
    sweep_sols(end+1, :) = {r, sol}; 
end

sol_prev = findSurfaceBVP2(r_init, A, n_mesh);
r = r_init;
while r > r_min + 1e-12
    target = max(r - r_step, r_min);

    bc = @(ya,yb,p) bcfun(ya, yb, p, target, A);
    try
        sol = bvp4c(@odefun, bc, sol_prev, opts);
        success = true;
    catch
        success = false;
    end

    if ~success
        step = (r - target)/2;
        ok = false;
        for k = 1:6
            candidate = r - step;
            bc_cand = @(ya,yb,p) bcfun(ya, yb, p, candidate, A);
            try
                sol = bvp4c(@odefun, bc_cand, sol_prev, opts);
                ok = true; target = candidate;
                break
            catch
                step = step/2;
            end
        end
        if ~ok
            fprintf('stalled near r=%.4f (backward)\n', r);
            break
        end
    end

    sol_prev = sol;
    r = target;
    sweep_sols = [{r, sol}; sweep_sols];   % prepend to keep sorted by r
end

%% produce the energies from the sweep solutions
% number of sweep points
n = size(sweep_sols,1);

% preallocate
E_bend = zeros(n,1);
E_int  = zeros(n,1);
R      = zeros(n,1);

% dimensionless bending energy and rim radius of each sweep shape
for ii = 1:n
    sol = sweep_sols{ii,2};     % grab the solution struct

    % bending energy: kappa*pi*E_bend = (kappa/2) * int (2H)^2 dA
    E_bend(ii) = int_bend_energy(sol);

    % endpoint radius r(1) from sol.y(3,end)
    r_end = sol.y(3,end);

    % E_int and R are the same number: the edge energy is gamma*2*pi*R*L0,
    % so E_int is used for the edge term and R for the perimeter axis
    E_int(ii) = r_end;
    R(ii)     = r_end;
end



%% Figure 2F
gamma_meas = 200;                  % line tension (kBT/um)
kappa_meas = 1200;                 % bending modulus (kBT)
rhog       = 0.04;                 % membrane weight per area (kBT/um^3)
A_meas = [1839.84, 3255, 4790];    % membrane areas (um^2)
L0 = sqrt(A_meas);
 
% cap energies 
%{ 
rad = A_meas./(2*pi*sqrt(A_meas/pi-(R.*L0).^2));
alpha = acos(A_meas./(2*pi.*rad.^2)-1);
phi = atan(2*tan(alpha/2));

E_bend_circ = 2*kappa_meas./rad.^2.*A_meas;
E_edge = 2*pi*gamma_meas*sqrt(A_meas./pi - A_meas.^2./(4*pi^2.*rad.^2));
E_grav_circ = 0.04*A_meas.*(-rad.*sin(alpha).*sin(phi) - rad/2.*(cos(alpha)+1).*cos(phi));
E_sphere = 4*pi*kappa_meas;
%}

figure; hold on
box on

% colours from the published Fig. 2F, one row per area in A_meas (blue, green, red)
lineCol = [0 0 245; 90 173 81; 234 52 36]/255;      % curves and marker edges
fillCol = [220 227 242; 229 240 219; 247 208 207]/255; % light marker fills

Astar =  1839.84232164833; %Astar = pi Rstar^2, Rstar = 2(2 kappa + kappabar)/gamma
crit_rhog = nan(1, length(A_meas)); % initialize vector of critial rhog

% start sweep
for i = 1:length(A_meas)
A_measi = A_meas(i); %i-th area in sweep

 % calculate grav energy
s_range = linspace(0,1,100);

nRows = size(sweep_sols,1);
L0 = sqrt(A_measi);
clear E_grav
for k = 1:nRows

        sol = sweep_sols{k,2};

        % evaluate r(s) and z(s)
        Y = deval(sol, s_range);   % 6×100
        r = Y(3,:);
        z = Y(4,:);

        rdata = [z z];
        zdata = [-r r];


        susp_pt = length(rdata);

        [x_rot, y_rot, phi, xcm_rot, ycm_rot] =  hang_curve([z z], [-r r], susp_pt);
       
        
        hvec(k) = y_rot(susp_pt)-ycm_rot;   % depth of the centre of mass below the hanging rim point (units of L0)
        h = hvec(k);
        E_grav(k,:) = A_measi.*h.*L0;
end

% total energy of the hanging axisymmetric shape (kBT):
%   kappa*pi*E_bend    = (kappa/2) * int (2H)^2 dA  (no kbar term: kbar = 0 in this model)
%   2*pi*gamma*L0*E_int = gamma * rim perimeter
%   rhog*E_grav         = rhog*A*h*L0 = weight * depth of the centre of mass
%                         below the hanging point (shapes computed without gravity)
E_axi_total = E_bend*kappa_meas*pi + 2*pi*gamma_meas.*sqrt(A_measi).*E_int - rhog.*E_grav; %total energy for axi model

x_ax = 2*pi*R*L0; %perimeter

% plot energy for perimeter > 2 um: below L ~ 0.5 um the sweep shapes depend
% on the bvp4c tolerance (see the sweep section), so tiny pores are left out
keep = x_ax > 2;
Lk = x_ax(keep);
Ek = E_axi_total(keep)/1e4;   % energy in units of 10^4 kBT
plot(Lk, Ek, '-.', 'LineWidth', 2, 'Color', lineCol(i,:))

% square = flat disk (largest L, last sweep point);
% circle = most closed shape plotted (smallest L > 2 um)
plot(Lk(end), Ek(end), 's', 'MarkerSize', 14, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', lineCol(i,:), 'MarkerFaceColor', fillCol(i,:))
plot(Lk(1), Ek(1), 'o', 'MarkerSize', 14, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', lineCol(i,:), 'MarkerFaceColor', fillCol(i,:))

% plot energy vs perimeter for caps
% plot((R*L0*2*pi), E_bend_circ + E_edge + E_grav_circ, 'LineStyle','-','LineWidth',2,'Color','k')

% reverse x-axis
set(gca,'XDir','reverse')

end

% axis ranges and labels as in the published Fig. 2F
xlim([-20 270])
ylim([2.3 4.5])
xlabel('$L\ (\mu\mathrm{m})$', 'Interpreter', 'latex')
ylabel('$E\ (10^4\,k_B T)$', 'Interpreter', 'latex')


%% Figure 3D
gamma_meas = 200;
kappa_meas = 1200;
A_meas = linspace(1000,38000,500);
L0 = sqrt(A_meas);


rhog = 0.04;   % 200 nm rods (kBT/um^3)
kbar = 20;     % saddle-splay modulus (kBT) 

% scales for the axes: x = A/A*, y = rhog*(2kappa + kbar)^2/gamma^3
K = 2*kappa_meas + kbar;
Rstar = 2*(2*kappa_meas + kbar)/gamma_meas;          % R* (um)
rhogfactor = (2*kappa_meas + kbar)^2/gamma_meas^3;   % y = rhog*rhogfactor
Astar = 4*pi*(2*kappa_meas + kbar)^2/gamma_meas^2;   % A* = pi*R*^2 (um^2)


figure; hold on
box on
grid on

% colours from the published Fig. 3D
purpleCol = [133 74 127]/255;   % axisymmetric curve where the barrier vanishes (crit_rhog)
goldCol   = [179 145 73]/255;   % spherical-cap curve where the barrier vanishes (Mathematica)
redCol    = [215 48 41]/255;    % 200 nm rods, rhog = 0.04
blueCol   = [16 7 231]/255;     % 200 nm rods, rhog = -0.04
greenCol  = [81 155 57]/255;    % 315 nm rods


% gold: spherical caps (from Mathematica). To second order in the cap
% curvature eps, the hanging-cap energy changes by
%   eps^2 * [ (2kappa+kbar)*A - gamma*A^(3/2)/(4*sqrt(pi)) + 3*rhog*A^(5/2)/(32*pi^(3/2)) ];
% the barrier vanishes where the bracket is zero, i.e. y = (2/3)*(x^-1 - 2*x^-3/2).
func = @(A) (2/3).*(A.^(-1)).^(5/2).*((-2).*A+A.^(3/2));
AA = linspace(0, 30, 1000);   % 1000 points (was 100) so the steep part is smooth
plot(AA, func(AA), 'linewidth',2, 'color', goldCol);

% black: a flat disk and a closed sphere (each hanging from one point) have
% equal energy,
%   gamma*2*sqrt(pi*A) - rhog*A*sqrt(A/pi) = 4*pi*(2kappa+kbar) - rhog*A*sqrt(A/(4*pi)),
% which in x = A/A* and y becomes 1 + y*x^(3/2) - x^(1/2) = 0.
func = @(Area , sigma) 1 + sigma .* (Area).^(3/2) - (Area).^(1/2);
f = fimplicit(func, [0, 50, -0.2,0.2]);
f.LineWidth=2;
f.Color = 'black';

% plot the rhog (2kappa + kappbar)^2/gamma^3 line for 315 nm rods and 200 nm rods
y1 = yline(rhog*rhogfactor);    % 200 nm rods
y2 = yline(0.063 * (2*5800 + 40)^2/500^3); % 315 nm rods: rhog = 0.063, kappa = 5800, kbar = 40, gamma = 500
y3 = yline(-rhog*rhogfactor);   % 200 nm rods with rhog = -0.04
y1.LineWidth = 2;
y1.LineStyle = ':';
y1.Color = redCol;
y1.Alpha = 1;       % yline is drawn 70% opaque by default
y2.LineWidth = 2;
y2.LineStyle = ':';
y2.Color = greenCol;
y2.Alpha = 1;
y3.LineWidth = 2;
y3.LineStyle = ':';
y3.Color = blueCol;
y3.Alpha = 1;

% set paper properties 
set(gcf, 'PaperSize', [5.5 4.5]);
set(gcf, 'PaperPosition', [0 0 5.5 4.5]);
set(gca, 'FontSize', 14)
set(gca, 'FontName', 'TimesNewRoman')


% sig figs and axis ticks 
set(gca,'XMinorTick','on','YMinorTick','on')
ytickformat('%.2f')
xticks(0:5:30)
xlim([0,12])        % published range (was [0,30])
ylim([-0.1, 0.16])  % published range (was [-0.2,0.2])

% axis labels 
xlabel('$A/A^*$', 'Interpreter','latex')
ylabel('$\rho g (2\kappa + \bar{\kappa})^2/ \gamma^3$', 'Interpreter','latex')


% purple: axisymmetric shapes. For each area, find the rhog at which
% dE/dR = 0 at the flat disk (the barrier vanishes), searching over rhogvec.
rhogvec = linspace(-0.24,0.1,2000);   % trial rhog values (kBT/um^3)

crit_rhog = nan(1, length(A_meas)); % initialize vector of critial rhog

% start sweep
for i = 1:length(A_meas)
A_measi = A_meas(i); %i-th area in sweep 

 % calculate grav energy
s_range = linspace(0,1,100);

nRows = size(sweep_sols,1);
L0 = sqrt(A_measi);
clear E_grav
for k = 1:nRows

        sol = sweep_sols{k,2};

        % evaluate r(s) and z(s)
        Y = deval(sol, s_range);   % 6×100
        r = Y(3,:);
        z = Y(4,:);

        rdata = [z z];
        zdata = [-r r];


        susp_pt = length(rdata);

        [x_rot, y_rot, phi, xcm_rot, ycm_rot] =  hang_curve([z z], [-r r], susp_pt);
       
        
        hvec(k) = y_rot(susp_pt)-ycm_rot;   % depth of the centre of mass below the hanging rim point (units of L0)
        h = hvec(k);
        E_grav(k,:) = A_measi.*h.*L0;
end

% same energy as in Figure 2F; E_grav is n x 1 and rhogvec is 1 x 2000, so
% E_axi_total is n x 2000 with one column per trial rhog
E_axi_total = E_bend*kappa_meas*pi + 2*pi*gamma_meas.*sqrt(A_measi).*E_int - rhogvec.*E_grav; %total energy for axi model




% eliminate the funny behaviour near R = 0 
[ind, ~] = find(R>0.1*R(end));

% differentiate each curve wrt R (each curve is for a diff rhog) 
derivatives = diff(E_axi_total(ind,:),1); 

% need to find the rhog that gives the stationary point dE/dR = 0
% As an approximation, just get all indices for dE/dR for each rhog that
% are positive. For these rhog the energy decreases as the flat disk
% (largest R) starts to curve, i.e. there is NO barrier; dE/dR < 0 at the
% disk means a barrier.
[ind1, ind2] = find(derivatives(end,:)>0);


% Take the rhog that gives the smallest positive dE/dR near the flat disk
% This is an approximation since we'd have to also find the rhog that
% gives dE/dR ~< 0. Then, the rhog that gives dE/dR = 0
% is between the rhog that gives dE/dR ~< 0 and the rhog that gives dE/dR
% ~> 0. But if we have a sufficiently large number of rhog's, we can just
% take the rhog that gives the smallest dE/dR ~> 0 as the approximation.  
[Mm, Ii] = min(derivatives(end,ind2)); 
crit = rhogvec(ind2(Ii));   % Ii counts within ind2, so map it back to rhogvec


if isempty(crit)
    continue
else 
    crit_rhog(i) = crit;
end

end

% shaded regions as in the published Fig. 3D:
%   gray   = open metastable: below the black curve and above the purple curve
%   purple = closed stable:   below the purple curve
Ahat  = A_meas/Astar;                        % 1 x 500
y_pur = crit_rhog*rhogfactor;                % purple curve (NaN where crit_rhog < min(rhogvec))
Ag    = linspace(min(Ahat), max(Ahat), 2000); % common grid for the two boundaries
y_blk = (sqrt(Ag) - 1)./Ag.^(3/2);           % black curve: 1 + y*A^(3/2) - A^(1/2) = 0
ok    = ~isnan(y_pur);
y_lo  = interp1(Ahat(ok), y_pur(ok), Ag);    % purple curve on the grid
y_lo(isnan(y_lo)) = -1;   % small A: purple curve lies below the plotted range, so extend to -1 (clipped by ylim)
y_hi  = max(y_blk, y_lo);
hG = fill([Ag, fliplr(Ag)], [y_hi, fliplr(y_lo)], [242 242 242]/255, 'EdgeColor', 'none');
hP = fill([Ag, fliplr(Ag)], [y_lo, -ones(size(Ag))], [232 224 240]/255, 'EdgeColor', 'none');
uistack([hG hP], 'bottom');   % shading behind all curves
set(gca, 'Layer', 'top')       % grid lines and ticks drawn on top of the shading

plot(Ahat, y_pur, 'LineWidth', 2, 'Color', purpleCol)   % purple: axisymmetric barrier-vanishing curve






function sol = findSurfaceBVP2(R1, A, n)
%FINDSURFACEBVP2  Solve the axisymmetric open-membrane shape via bvp4c
%
%   sol = findSurfaceBVP2(R1,A,n) finds the minimum-bending-energy shape
%   with total area A and rim radius R1 (lengths in units of L0 =
%   sqrt(area), so A = 1), without gravity, starting from n mesh points.
%   It solves dy/dx = f(x,y,p), x in [0,1], for y = [a; da; r; z; q; A]
%   with unknown parameters p = [L; sigma] (see odefun and bcfun).
%
%   Example:
%     sol = findSurfaceBVP2(0.2, 1.0, 20);
%     xFine = linspace(0,1,400);
%     yFine = deval(sol,xFine);

% parameter guess [L; sigma]
pinit   = [1.0; 0.0];

% constant‐vector guess function 
guessFcn = @(x)[ ...
    pi/2;    ... % a(x)
    1;       ... % a'(x)
    R1;      ... % r(x)
    1;       ... % z(x)
    0;       ... % q(x)
    A/2        % A(x)
    ];

% initial mesh
xmesh = linspace(0,1,n);

% build solver‐init struct 
solinit = bvpinit(xmesh, guessFcn, pinit);
opts = bvpset('RelTol',1e-6, 'AbsTol',1e-6, 'Stats','on');

% solve 
sol = bvp4c(@odefun, @(ya,yb,p)bcfun(ya,yb,p,R1,A), solinit, opts);

end

function dydx = odefun(~, y, p)
% Axisymmetric shape equations with kappa = 1, no pressure and no gravity
% (Lagrange-multiplier form, cf. Juelicher & Seifert 1994). x in [0,1] is
% arclength divided by the total length L, so d/ds = (1/L) d/dx, and
% a = psi is the tangent angle (r_s = cos psi, z_s = sin psi):
%   psi_ss = -cos(psi)*psi_s/r + sin(psi)*cos(psi)/r^2 + q*sin(psi)/r
%   r_s = cos(psi),   z_s = sin(psi),   A_s = 2*pi*r
%   q_s = sigma + (psi_s^2 - sin(psi)^2/r^2)/2
% q is the Lagrange multiplier that enforces r_s = cos(psi), and sigma is
% the tension (Lagrange multiplier for the total area).
% y = [a; da; r; z; q; A], with da = d(psi)/dx
L     = p(1);
sigma = p(2);

a  = y(1,:);
da = y(2,:);
r  = y(3,:);
z  = y(4,:);  % not used below
q  = y(5,:);
A  = y(6,:);  % not used below


dydx = [ ...
    da; ...
    -L.*cos(a)./r .* da + L.^2.*(sin(2.*a)./(2.*r.^2) + q.*sin(a)./r); ...
    L.*cos(a); ...
    L.*sin(a); ...
    L.*sigma + (da.^2./L - L.*sin(a).^2./r.^2)/2; ...
    2*pi*L.*r ...
    ];



end

function res = bcfun(y0, y1, p, R1_, A_)
L = p(1);


% left endpoint (x=0)
a0  = y0(1);  da0 = y0(2);
r0  = y0(3);  z0  = y0(4);
q0  = y0(5);  A0  = y0(6);

% right endpoint (x=1)
a1  = y1(1);  da1 = y1(2);
r1  = y1(3);   
A1  = y1(6);

% 8 BCs for 6 unknowns + 2 params. The rim radius is prescribed, so the
% line tension does not enter the shape, only the energy.
res = [ ...
    z0;                    % z(0)=0
    r0 - 0.0001;           % r(0)=1e-4: pole, just off the axis to avoid 1/r
    r1 - R1_;              % r(1)=R1: rim radius
    A0;                    % A(0)=0
    a0;                    % a(0)=0: smooth pole
    A1 - A_;               % A(1)=A: total area
    r1*da1 + L*sin(a1);    % r(1)a'(1)+L sin(a(1))=0, i.e. 2H = 0 at the free rim (no bending moment, kbar = 0)
    q0                     % q(0)=0: no force at the pole
    ];
end

function E_bend = int_bend_energy(sol)
% Dimensionless bending energy E_bend = int_0^1 (L sin(psi) + r psi')^2/(L r) dx.
% Since 2H = psi_s + sin(psi)/r = (r psi' + L sin(psi))/(L r) and
% dA = 2*pi*r*L dx, kappa*pi*E_bend = (kappa/2) * int (2H)^2 dA.
% It is independent of the length scale L0.

L = sol.parameters(1);

s_integrand = @(s) localIntegrand(sol,L,s);

E_bend = integral(s_integrand, 0, 1,  'RelTol',1e-6, 'AbsTol',1e-9);

end

function val = localIntegrand(sol, L, s)
Y = deval(sol, s);          % 6 x numel(s)

a  = Y(1,:);
da = Y(2,:);
r  = Y(3,:);

val = (L.*sin(a) + r.*da).^2 ./ (L .* r);
end


function [x_rot,y_rot,phi,xcm_rot,ycm_rot,dx_res] = hang_curve(x,y,k)
% Rotate a 2D curve (surface of revolution) given its points x,y and pivot
% index k 
% Returns the rotated points and the rotated COM 

% make into column vector 
x = x(:); y = y(:); 

% compute COM 
[xcm, ycm] = com_surface_yaxis(x,y); 

% choose pivot point
xs = x(k); ys = y(k);

% create ray from pivot point to COM point
dx = xcm - xs;  dy = ycm - ys;
% this angle is the angle the ray makes wrt the + y-axis where the left is
% negative and the right is positive in (-pi,pi]. We want the angle between
% the -y axis and the ray, so that we can rotate the ray and make it point
% towards -pi/2. To do so, we first modulo 2*pi and then subtract pi. This
% gives not only the angle between the ray and -y axis, it also gives the
% correct sign for our rotation: For a ray pointing to the left downward
% direction wrt the + y-axis, with angle -2pi/3, we want to rotate that ray
% CCW. The result of the below gives pi/3, which will rotate the ray CCW to
% -pi/2 with the rotation matrix in rot_about. 
phi = atan2(dx,dy);
phi = mod(phi,2*pi) - pi; 

% rotate curve
[x_rot,y_rot]     = rot_about(x,y,xs,ys,phi);
% rotate COM 
[xcm_rot,ycm_rot] = rot_about(xcm,ycm,xs,ys,phi);

% sanity check: is the rotated COM below the pivot? 
dx_res = xcm_rot - xs;                   
end

function [xr,yr] = rot_about(x,y,xc,yc,phi)
% this code performs the rotation

% create rays for every surface point from the pivot point xc, yc 
u = x - xc; v = y - yc; c = cos(phi); s = sin(phi);

% rotate those rays with the angle phi (the sign in phi gives the CCW or CW
% direction) and then reshift the rotated rays back to the pivot pt xc,yc
xr =  c.*u - s.*v + xc;
yr =  s.*u + c.*v + yc;
end


function [x_cm, y_cm] = com_surface_yaxis(x, y)

% this code computes the COM for a surface of revolution
% A = 2*pi*int y ds ~= 2 pi sum (y * ds)
% Note: with x = [z z], y = [-r r] the segment joining the two halves (rim
% back to pole) is also counted; this changes h by < 1% (checked against the
% exact centre of mass of the BVP solution).

x = x(:); y = y(:);                                  % make into column vector 

dx      = diff(x); dy = diff(y); ds = hypot(dx,dy);  % differentiate, and compute ds = sqrt(dx^2 + dy^2)
xm      = 0.5*(x(1:end-1)+x(2:end));                 % compute midpoints in x
ym      = 0.5*(y(1:end-1)+y(2:end));                 % compute midpoints in y
w       = 2*pi*abs(ym).*ds;          
y_cm    = 0;                                         % by symmetry
x_cm    = sum(w .* xm) / sum(w);
end



function [H, ks, kphi, s, psi] = axisym_mean_curvature(s,psi,psip)
% AXISYM_MEAN_CURVATURE  Mean curvature H = (ks + kphi)/2 of an axisymmetric
% profile, from arclength s (um), tangent angle psi and psip = dpsi/ds (1/um).
%   ks   = dpsi/ds       (meridional curvature)
%   kphi = sin(psi)/r    (azimuthal curvature), with r rebuilt here as
%                        int cos(psi) ds starting from r = 0
%
% Outputs
%   H    : mean curvature (size N)
%   ks   : meridional principal curvature (along the profile)
%   kphi : azimuthal principal curvature (around the hoop)
%   s    : cumulative arclength parameter (starts at 0)
%   psi  : tangent angle, where r_s = cos(psi), z_s = sin(psi)
%



r0 = 0; z0 = 0;
% Reconstruct profile from psi(s)
r = r0 + cumtrapz(s, cos(psi));
z = z0 + cumtrapz(s, sin(psi));


ks = psip;

tol_axis = 5 * eps(max(r));
kphi = sin(psi) ./ max(r, tol_axis);
near_axis = abs(r) < tol_axis;
if any(near_axis)
    kphi(near_axis) = ks(near_axis);
end



% Mean curvature
H = 0.5 * (ks + kphi);
H(1) = ks(1);


end

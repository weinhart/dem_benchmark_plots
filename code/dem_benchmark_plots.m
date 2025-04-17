%% DEM Benchmark – Post‑processing & Plotting Script
% -------------------------------------------------------------------------
% This script reproduces the key figures from the benchmark paper
% "Comparing open‑source DEM frameworks for simulations of common bulk
% processes" (Computer Physics Communications 296 (2024) 109066).
% Figures in the paper ↔ files written by this script:
%   • Fig. 6   – Velocity histograms in the silo        ➜ ../output/figure6.png
%   • Fig. 9   – Centre of mass trajectories (drum)    ➜ ../output/figure9.png
%   • Fig. 11  – Penetration depth of steel impactor   ➜ ../output/figure11.png
% 
% After cloning the repository simply run:
%     >> dem_benchmark_plots
% in the MATLAB console.  All figures will be regenerated in ../output/.
% 
% Estimated runtimes on a laptop: < 1 min.
% -------------------------------------------------------------------------
% Author : T Weinhart – University of Twente
% Created: 2023‑04‑01 – last edit 2025‑04‑17
% -------------------------------------------------------------------------

%% 0. Clear workspace 
fprintf('[%s] Clean start – clearing variables\n', datestr(now,'HH:MM:SS'))
clearvars           % remove variables only
close all           % close open figures
clc                 % clear command window

%% 1. Set colour palette & solver names ----------------------------------
fprintf('Assigning colours and solver display names...\n');
% Obtain a visually distinct 9‑colour palette (requires linspecer.m)
c = linspecer(9);

% Assign one colour per solver (indices match palette used in the paper)
MercuryDPM.c = c(6,:);
MUSEN.c      = c(8,:);
MFiX.c       = c(7,:);
Kratos.c     = c(4,:);
LIGGGHTS.c   = c(5,:);
Kratos.c     = c(4,:);   % duplicate assignment kept to preserve original code
Blaze.c      = c(1,:);
Granoo.c     = c(3,:);
Yade.c       = c(9,:);

% Human‑readable solver names (appear in legends)
MercuryDPM.name = "MercuryDPM";
MUSEN.name      = "MUSEN";
MFiX.name       = "MFiX";
Kratos.name     = "Kratos";
LIGGGHTS.name   = "LIGGGHTS-public";
Blaze.name      = "Blaze-DEM";
Granoo.name     = "GranOO";
Yade.name       = "Yade";

%% 2. Load drum COM data --------------------------------------------------
fprintf('Loading drum centre‑of‑mass data...\n');
drumRaw = xlsread('../data/ResultsRev1DrumCOM.xlsx',1,'B:ZZ');

% Five‑column blocks per solver: [time x z x_flipped z_flipped]
MercuryDPM.drum = drumRaw(:,[1 4 5 2 3]);
MUSEN.drum      = drumRaw(:,6:10);
MFiX.drum       = drumRaw(:,11:15);
Kratos.drum     = drumRaw(:,16:20);
LIGGGHTS.drum   = drumRaw(:,21:25);
Blaze.drum      = drumRaw(:,26:30);
Granoo.drum     = drumRaw(:,31:35);

% Flip X/Z axes for Granoo (mirrored coordinate system)
fprintf('Flipping X/Z axes for Granoo data...\n');
Granoo.drum(:,4) = -Granoo.drum(:,4);
Granoo.drum(:,2) = -Granoo.drum(:,2);

% Yade block plus axis flip
Yade.drum       = drumRaw(:,36:40);
Yade.drum(:,4)  = -Yade.drum(:,4);
Yade.drum(:,2)  = -Yade.drum(:,2);

%% 3. Load penetration test data -----------------------------------------
fprintf('Loading steel‑ball penetration data...\n');
penRaw = xlsread('../data/ResultsRev1Penetration25K.xlsx',1,'B:ZZ');

MercuryDPM.pen = penRaw(:,1:11); % read 11 columns (time + data)
MercuryDPM.pen(~isfinite(MercuryDPM.pen(:,1)),:) = [];% remove empty rows
MUSEN.pen      = penRaw(:,11+(1:11));
MFiX.pen       = penRaw(:,22+(1:11));
Kratos.pen     = penRaw(:,33+(1:11));
LIGGGHTS.pen   = penRaw(1:20:end,44+(1:10)); % down‑sample rows
Blaze.pen      = penRaw(:,55+(1:11));
Granoo.pen     = penRaw(:,66+(1:11));
Yade.pen       = penRaw(:,77+(1:11));

%% 4. Load silo‑flow data (large orifice, material M1) -------------------
fprintf('Loading silo discharge data (large M1)...\n');
largeM1 = xlsread('../data/ResultsRev1SiloLargeM1.xlsx',1,'B:ZZ');
MercuryDPM.largeM1 = extractSiloProperties(largeM1(:,1:7));
MUSEN.largeM1      = extractSiloProperties(largeM1(:,[8 9 10 12 13 14 11]));
MFiX.largeM1       = extractSiloProperties(largeM1(:,7+[8 9 10 12 13 14 11]));
Kratos.largeM1     = extractSiloProperties(largeM1(:,14+[8 9 10 12 13 14 11]));
LIGGGHTS.largeM1   = extractSiloProperties(largeM1(:,21+[8 9 10 12 13 14 11]));
Blaze.largeM1      = extractSiloProperties(largeM1(:,28+[9:14 8]));
Granoo.largeM1     = extractSiloProperties(largeM1(:,35+[8:14]));
Yade.largeM1       = extractSiloProperties(largeM1(:,42+[8:14]));

%% 5. Drum COM plots ------------------------------------------------------
fprintf('Plotting drum COM trajectories...\n');
figure(2); clf;
set(gcf,'Position',[0 0 740 280]);         % same figure size as original
t = tiledlayout(1,2);                      % create two side‑by‑side axes

% --- Sub‑plot 1: flipped coordinate frame -------------------------------
nexttile; hold on;
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN, Yade}
    data = data_{1};
    ix   = data.drum(:,1) < 4.91;          % restrict to first drum revolution
    plot(data.drum(ix,4), data.drum(ix,5), '.', 'DisplayName', data.name, 'Color', data.c);
end
hold off;
axis equal;
axis([-0.0223 0.0701 -0.0782 -0.0168]);
xlabel('$x\,[m]$','Interpreter','latex');
ylabel('$z\,[m]$','Interpreter','latex');
set(gca,'Box','on');
h = legend('show');
set(h,'Box','off','Location','NorthWest');
h.Position = h.Position + [-0.02 0 0 0];    % manual tweak due to paper layout

% --- Sub‑plot 2: original frame -----------------------------------------
nexttile; hold on;
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN, Yade}
    data = data_{1}; ix = data.drum(:,1) < 4.91;
    plot(data.drum(ix,2), data.drum(ix,3), '.', 'DisplayName', data.name, 'Color', data.c);
end
hold off;
axis equal; axis tight;
axis([-0.0223 0.0701 -0.0782 -0.0168]);
xlabel('$x\,[m]$','Interpreter','latex');
ylabel('$z\,[m]$','Interpreter','latex');
exportgraphics(gcf,'../output/figure9.png','Resolution',500);

%% 6. Penetration depth plots --------------------------------------------
fprintf('Plotting penetration depth curves...\n');
figure(5); clf;
set(gcf,'Position',[0 0 500 250]);

t = MercuryDPM.pen(:,1);  % common time vector
m = []; s = [];            % pre‑allocate mean & std matrices

% --- Compute mean/std envelope across solvers (excluding Yade) ----------
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN}
    data = data_{1};
    data.pen = data.pen(isfinite(data.pen(:,1)),:);          % remove NaN rows
    m(:,end+1) = mean(interp1(data.pen(:,1), data.pen(:,2:end), t, [], 'extrap'),2);
    s(:,end+1) = std( interp1(data.pen(:,1), data.pen(:,2:end), t, [], 'extrap'),[],2);
end
m = mean(m,2) - 0.06;                      % create global vertical offset (change of variable)
s = sqrt(mean(s.^2,2));                    % std‑dev
up = [t, m+s]; low = [t, m-s]; comb = [up; flipud(low)];
fill(comb(:,1), comb(:,2), 0.95*[1 1 1], 'EdgeColor', [0.95 0.95 0.95], 'DisplayName','mean\pmstd');
hold on;

% --- Overlay each solver trajectory -------------------------------------
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN, Yade}
    data = data_{1};
    plot(data.pen(:,1), mean(data.pen(:,2:end),2)-0.06, '.-', 'DisplayName', data.name, 'Color', data.c);
end
hold off;
xlabel('Time (s)','Interpreter','latex');
ylabel('Vertical displacement (m)','Interpreter','latex');
legend('show');
set(legend,'location','north','Box','off');
set(gca,'Box','on'); ylim([-0.15 0]);
exportgraphics(gcf,'../output/figure11.png','Resolution',500);

%% 7. Silo velocity histograms -------------------------------------------
fprintf('Plotting silo velocity distributions...\n');
figure(6); clf;
set(gcf,'Position',[0 0 500 500]);

subplot(3,1,1); hold on;   % --- Vr histogram ---
N = 100; v = MercuryDPM.largeM1.Vr; bins = linspace(quantile(v,0.05), quantile(v,0.99), N);
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN, Yade}
    data = data_{1}.largeM1;
    histogram(data.Vr, bins, 'DisplayName', data_{1}.name, 'DisplayStyle','stairs', 'EdgeColor', data_{1}.c);
end
hold off; xlabel('$v_r\,[m/s]$','Interpreter','latex'); ylabel('Frequency','Interpreter','latex'); axis tight;
legend('show'); set(legend,'location','northwest','Box','off'); set(gca,'Box','on');

subplot(3,1,2); hold on;   % --- Oz histogram ---
v = MercuryDPM.largeM1.Oz; bins = linspace(-0.5,0.5,N);
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN, Yade}
    data = data_{1}.largeM1; histogram(data.Oz, bins, 'DisplayStyle','stairs', 'EdgeColor', data_{1}.c);
end
hold off; xlabel('$v_\theta\,[rad/s]$','Interpreter','latex'); ylabel('Frequency','Interpreter','latex'); axis tight; set(gca,'Box','on');

subplot(3,1,3); hold on;   % --- Vz histogram ---
v = MercuryDPM.largeM1.Vz; bins = linspace(quantile(v,0.05), 0, N);
for data_ = {Blaze, Granoo, Kratos, LIGGGHTS, MercuryDPM, MFiX, MUSEN, Yade}
    data = data_{1}.largeM1;
    fprintf('%15s %5d particles plotted\n', data_{1}.name, size(data.X,1));
    histogram(data.Vz, bins, 'DisplayStyle','stairs', 'EdgeColor', data_{1}.c, 'DisplayName', data_{1}.name);
end
hold off; xlabel('$v_z\,[m/s]$', 'Interpreter','latex'); ylabel('Frequency','Interpreter','latex'); axis tight; set(gca,'Box','on');
exportgraphics(gcf,'../output/figure6.png','Resolution',500);

%% 8. Helper function -----------------------------------------------------
function data = extractSiloProperties(raw)
% extractSiloProperties – cleans raw silo data block & computes derived vars
% INPUT : raw – numeric matrix with columns [X Y Z Vx Vy Vz dummy]
% OUTPUT: struct with Cartesian velocities + radial speed & angular vel.

raw(~isfinite(raw(:,1)),:) = [];          % remove NaN rows (empty lines)
raw(raw(:,3) < -0.06,:)   = [];           % discard particles below orifice

% Direct assignments ------------------------------------------------------
[X,Y,Z, Vx,Vy,Vz] = deal(raw(:,1), raw(:,2), raw(:,3), raw(:,4), raw(:,5), raw(:,6));

% Derived kinematics ------------------------------------------------------
Oz = (X .* Vy - Y .* Vx) ./ (X.^2 + Y.^2);               % Angular velocity
Vr = (Vx .* X  + Vy .* Y) ./ sqrt(X.^2 + Y.^2);          % Radial velocity

% Pack into struct --------------------------------------------------------
data = struct('X',X,'Y',Y,'Z',Z,'Vx',Vx,'Vy',Vy,'Vz',Vz,'Oz',Oz,'Vr',Vr);
end

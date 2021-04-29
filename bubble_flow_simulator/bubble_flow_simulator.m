clear all
close all
clc
%% Path and Folders Management
addpath(genpath('..\..\hatimb-particle_flow_simulator\'));
root = dir('..\..\');
root_dir = root(1).folder;
save_dir = '\hatimb-particle_flow_simulator_DATA';
mkdir(root_dir,save_dir);
save_path = [root_dir save_dir '\'];
%% Temp folder
mkdir(root_dir,[save_dir '\temp']);
%% Visualize
display = 0; % 0: No display, 1: Minimal display, 2: All displays
%% Essential Variables
disp('Running...')
samp_freq = 1000; %(Hz)
n_bubbles = 5000; % Number of bubbles trajectories generated
bb_per_paquet = 5000; % n_trajectories per paquet (for storage)
n_bubbles_steady_state = 100; % 1/5th is taken to avoid shortage of bubbles
 % 1/5th is taken to avoid shortage of bubbles
t_steady_state = 1; % Desired simulation time (s)
bubble_size = 2; % Bubble diameter (um)
pulsatility = 1; % 1 = Yes | 0 = No
file_name = 'test';
%% Loading
name = 'tree5'; % Name of the .swc graph model
filename = [name '.swc'];
g = importdata(filename);
%% Variables
target = g(:,1);    % Nodes IDs
source = g(:,7);    % Parent node
pos = g(:,3:5);     % Positions [x,y,z]
r = g(:,6);         % Nodes Radii
r_norm = r./max(r); % Normalized Radii [0,1]
r_inverse = 1./r;
r_inverse_norm = 1./r_norm;
%% Radii hist
if display
    figure(77);clf
    hist(r);xlabel('radius');ylabel('N');
    [counts,centers] = hist(r);
end
%% Viewing graph nodes
if display == 1
    figure(1);clf
    h = scatter3(pos(:,1),pos(:,2),pos(:,3),1,[1 1 1]); % Shortest path nodes);
    alpha = 0.2;
    set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
    darkBackground(gcf)   
    axis equal;
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('z (\mum)')
end
drawnow
%% Positions scaling % Use only for anisotropic dataset
% dim_1 = 2; % (um)
% dim_2 = 2; % (um)
% dim_3 = 50; % (um)
% pos(:,1) = pos(:,1) * dim_1;
% pos(:,2) = pos(:,2) * dim_2;
% pos(:,3) = pos(:,3) * dim_3;
% r = r .* (dim_1+dim_2)/2;
%% Graph
s = source+2;
t = target+2;
C = setxor(s,1:length(s)); % Finding missing parents
end_nodes = [C-1;s(end)];  % Adding last extremity
clear C
%% Finding bifurcations
[uniqueA, i, j] = unique(s,'first');        % Finding unique nodes
tmp = find(not(ismember(1:numel(s),i)));    % idx of duplicates
tmp2 = s(tmp);                              % duplicates
[biff_nodes, ii, jj] = unique(tmp2,'first');% Keeping only first occurence
biff_nodes = biff_nodes-1; % it is the previous node
if display
    figure(72);
    clf
    %scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled');
    plot3(pos(:,1),pos(:,2),pos(:,3),'.k');
    hold on
    plot3(pos(biff_nodes,1),pos(biff_nodes,2),pos(biff_nodes,3),'og') % Starting flow node
    title('Bifurcations and endnodes');
    plot3(pos(end_nodes,1),pos(end_nodes,2),pos(end_nodes,3),'or') % Starting flow node
end
%% Directed graph visualisation
if display == 1
    disp('Directed graph visualisation...')
    start = 1;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
    finish = start + randi([1 size(s,1)-start-1],1,1); % random integer from [start:#target_nodes-start]
    DG = digraph(s,t); % Directed graph generation
    [SP, D] = shortestpathtree(DG,start,finish); % Shortest path
    edges = table2array(SP.Edges); % Conversion
    nodes = [edges(:,1);edges(end,2)];
    trajectory = pos(nodes,:);
    %Plot 1 trajectory
    figure(2);clf;f = plot(SP); % Plot the shortest path graph variable
    f.XData = [pos(1,1);pos(:,1)];
    f.YData = [pos(1,2);pos(:,2)];
    f.ZData = [pos(1,3);pos(:,3)];
    f.EdgeColor = 'r';
    f.LineWidth = 3;
    f.ArrowSize = 10;
    hold on
    plot3(pos(biff_nodes,1),pos(biff_nodes,2),pos(biff_nodes,3),'og') % Starting flow node
    plot3(pos(end_nodes,1),pos(end_nodes,2),pos(end_nodes,3),'ok') % Starting flow node
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    zlabel('z','FontSize',20);
    view(20,30)
    set(gcf,'color','w');
end
%% Volume calculation
total_vessels_length = 0;
vectors_volume_calc = pos(t(1:length(t)-1),:)-pos(s(1:length(t)-1),:); % Parallel vectors between all nodes
distances_volume_calc = sqrt(sum(diff(vectors_volume_calc,[],1).^2,2)); % Euclidian norms calculation
length_volume_calc = sum(distances_volume_calc);% um
mean_radius = mean(r); % um
total_vessel_volume = length_volume_calc*pi()*mean_radius.^2; % um^3 % This is the cumulative trajectories volume! Not network volume
% disp(['Total vasculature volume = ' num2str(total_vessel_volume*1E-9) ' mm^3']);
whole_volume = (max(pos(:,1))-min(pos(:,1)))*(max(pos(:,2))-min(pos(:,2)))*(max(pos(:,3))-min(pos(:,3)));
% disp(['Vessel volume ratio = ' num2str(total_vessel_volume/whole_volume*100) '%']);
vessel_volume_mL = total_vessel_volume*1E-9/1000;
C_MB = 2E5; % Concentration of MicroBubbles in the bloodstream, Hingot et al, 2019
n_bubbles_in_network = round(vessel_volume_mL*C_MB); % Number of microbubbles in the studied vessels
% disp(['Number of MB in vessels : ' num2str(n_bubbles_in_network)]);
%% Extremity points visualization (Computationally heavy)
if display == 3
    disp('Extremity points visualization')
    figure(3);clf;
    plot3(pos(:,1),pos(:,2),pos(:,3),'.k');
    hold on
    plot3(pos(end_nodes,1),pos(end_nodes,2),pos(end_nodes,3),'or');
    plot3(pos(biff_nodes,1),pos(biff_nodes,2),pos(biff_nodes,3),'og') % Starting flow node
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    zlabel('z','FontSize',20);
    legend('Nodes','Endnodes','Bifurcations','location','best')
    clear DG_display
end
%% Vincent Hingot stats
%%% N
log_d = linspace(0,5,1000);
log_N = 3.7*log_d -8.2;
d = exp(log_d);
N = exp(log_N);
log_v = 1.9*log_d -6;
v = exp(log_v);
d_sample_log = log(2*r);
v_sample_log = 1.9*(d_sample_log) -6;
v_sample = exp(v_sample_log);
v_sample_um = v_sample*1000;
d_sample = 2*r;
N_sample_log = 3.7*d_sample_log -8.2;
N_sample = exp(N_sample_log);
%%% Affichage
if display == 2 
    figure(4);clf
    subplot(2,3,1);plot(log_d,log_N,'LineWidth',2);
    hold on; plot(d_sample_log,N_sample_log,'.');
    grid on;title('Dependency of the bubble rate with vessel’s diameter');xlabel('log(d)');ylabel('log(N)');
    axis([0 5 0 8]);
    subplot(2,3,2);plot(d,N,'LineWidth',2);
    hold on; plot(d_sample,N_sample,'.');
    grid on;title('Dependency of the bubble rate with vessel’s diameter');xlabel('d (um)');ylabel('N');
    subplot(2,3,3);plot(d_sample,N_sample,'.');
    grid on;title('Sample dependency of the bubble rate with vessel’s diameter');xlabel('d (um)');ylabel('N');
    %%% v
    subplot(2,3,4);plot(log_d,log_v,'LineWidth',2);
    hold on; plot(d_sample_log,v_sample_log,'.');
    legend('Hingot','Sample data','Location','Best');
    grid on;title('Dependency of maximum velocity with vessel’s diameter');xlabel('log(d)');ylabel('log(v)');
    axis([0 5 -5 4]);
    subplot(2,3,5);plot(d,v,'LineWidth',2);
    hold on; plot(d_sample,v_sample,'.');
    legend('Hingot','Sample data','Location','Best');
    grid on;title('Dependency of maximum velocity with vessel’s diameter');xlabel('d (um)');ylabel('v (mm/s)');
    subplot(2,3,6);plot(d_sample,v_sample,'.');
    grid on;title('Dependency of maximum velocity with vessel’s diameter');xlabel('d (um)');ylabel('v (mm/s)');
    %figure;plot(d,v,'LineWidth',2);xlabel('d (um)');ylabel('v');title('Sample velocities');
    %axis([15 45 0 3.5]);
end
%% Poiseuille distribution
x = linspace(-1,1,1000);
v_poiseuille = 1-x.^2;
v_poiseuille_squared = v_poiseuille.^2;
% figure;plot(x,v_poiseuille);hold on;plot(x,v_poiseuille_squared);
% %legend('Poiseuille velocity','Velocity squared','Location','Best');
% xlabel('Normalized diameter')
% ylabel('Poiseuille value (P)')
% title('Poiseuille Distribution');
%% Pulsatility function
BPM = 300;
freq = BPM/60;
period = 1/freq;
dt = 1/samp_freq;
t_f = 200;
if t_steady_state > t_f
    t_f = t_steady_state;
end
x = 0:dt:t_f-dt;
ecg_raw = ecg(BPM,dt,t_f);
ecg_filtered = ecg_raw-min(ecg_raw);
ecg_filtered2 = ecg_filtered./max(ecg_filtered);
ecg_filtered3 = ecg_filtered2+0.5;
ecg_normalized = ecg_filtered3;
% figure;plot(x,ecg_normalized);xlabel('time(s)');ylabel('Normalized Amplitude');
% figure(4)
% clf
% %plot(ecg_raw)
% hold on 
% % plot(x ,ecg_filtered)
% % plot(x ,ecg_filtered2);
% plot(x ,ecg_filtered3,'LineWidth',1.5);
% xlabel('Time (s)');
% ylabel('Multiplication Factor');
% title([num2str(BPM) ,' BPM']);
% % legend('ecg filtered','ecg filtered2','ecg filtered3','location','best');
% set(gca,'FontSize',14)
% grid on
clear ecg_filtered3 ecg_filtered2 ecg_filtered ecg_raw
%% Trajectories statistics
% From start to endnodes
disp('Computing trajectories statistics...');
DG = digraph(s,t,r_inverse); % Directed graph generation
end_nodes_biff = [end_nodes;biff_nodes];%s(2:end);%[end_nodes;biff_nodes]; % Merge endnodes and bifurcation nodes
d_TRAJECTORIES = zeros(1,numel(end_nodes_biff));
mean_RADII = zeros(1,numel(end_nodes_biff));
median_RADII = zeros(1,numel(end_nodes_biff));
min_RADII = zeros(1,numel(end_nodes_biff));
max_RADII = zeros(1,numel(end_nodes_biff));
tic
for idx = 1:numel(end_nodes_biff)
    start = 1;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
    [SP, ~] = shortestpathtree(DG,start,end_nodes_biff(idx)); % Shortest path
    edges = table2array(SP.Edges);
    nodes = [edges(:,1);edges(end,2)]; % It is the previous node!
    trajectory = pos(nodes,:); % Nodes positions attribution
    d_TRAJECTORIES(idx) = sum(sqrt(sum(diff(trajectory,[],1).^2,2))); % Total length
    mean_RADII(idx) = mean(r(nodes));
    median_RADII(idx) = median(r(nodes));
    min_RADII(idx) = min(r(nodes));
    max_RADII(idx) = max(r(nodes));
    if idx==1
       time_estimate = toc*numel(end_nodes_biff);
       fprintf('It should take ~%1.0f min.\n',time_estimate/60);
    end
end
[mean_RADII_sorted,Idx_mean] = sort(mean_RADII,'descend');
[median_RADII_sorted,Idx_median] = sort(median_RADII,'descend');
[min_RADII_sorted,Idx_min] = sort(min_RADII,'descend');
[max_RADII_sorted,Idx_max] = sort(max_RADII,'descend');
d_TRAJECTORIES_norm = min(d_TRAJECTORIES)./d_TRAJECTORIES;
if display == 3
    figure;clf
    plot(d_TRAJECTORIES,'.');title('Length');ylabel('Trajectry length (\mum)')
    figure;clf
    subplot(1,4,1);
    plot(mean_RADII_sorted,'.');title('Radius - Mean');ylabel('Mean trajectory radius (\mum)');
    subplot(1,4,2);
    plot(median_RADII_sorted,'.');title('Radius - Mean');ylabel('Mean trajectory radius (\mum)');
    subplot(1,4,3);
    plot(min_RADII_sorted,'.');title('Radius - MIN');ylabel('Min trajectory radius (\mum)');
    subplot(1,4,4);
    plot(max_RADII_sorted,'.');title('Radius - MAX');ylabel('Max trajectory radius (\mum)');
end
toc
%% Trajectories selection probability
disp('Computing trajectories selection probability...');
radii = min_RADII_sorted; % rounding to units
% end_nodes_sorted = end_nodes(Idx_min);
end_nodes_biff_sorted = end_nodes_biff(Idx_min);
radii_rounded = round(min_RADII_sorted);
radii_unique = unique(radii_rounded);
% radii_unique_continuous = max(radii_unique):-1:min(radii_unique);
n_radii = numel(radii_unique);% number of differrent radii
slope = 9; % %%% Tweaking to a higher slope to get close to Hingot's
intercept = 0;
N_traject_log = slope*(log(radii*2))+intercept; %%% TEMPORARY FOR TESTING
% N_traject_continuous_log = 3.7*(log(radii_unique_continuous*2)) -8.2;
N_traject = exp(N_traject_log);
(log(N_traject(1))-log(N_traject(2)))/(log(radii(1)*2)-log(radii(2)*2));
N_traject(1);
N_traject(2);
% N_traject_continuous = exp(N_traject_continuous_log);
% N_traject_continuous_norm = N_traject_continuous/sum(N_traject_continuous);
% Let's normalize the N so that sum(N) = 1
% Let's account for that
radii_count = histc(radii_rounded, radii_unique); % this willgive the number of occurences of each unique element
radii_count = fliplr(radii_count);
radii_unique = sort(radii_unique,'descend');
for i = 1:numel(radii_unique)
    start = sum(radii_count(1:i-1)) + 1;
    finish = sum(radii_count(1:i));
    N_traject_norm(start:finish) = N_traject(start:finish)/radii_count(i);
end
% N_traject_norm = N_traject_norm.*(d_TRAJECTORIES_norm.^3.5); % compensate probability with length
N_traject_norm = N_traject_norm/sum(N_traject_norm); % normalize for pdf according to trajectory length
%% Simuation
disp('Starting simulation...');
padding_bubble = bubble_size/2; % To account for the fact that the bubbles are not infinitesimal points
% bubbles = cell(1,1);%cell(n_bubbles,1); % Initialization
tot_toc = 0; % For displaying progress to the user
min_poiseuille = 0.2; % Minimum Poiseuille value (a value of 0 causes an infinite computation time since the bubble doesn't move)
min_length = 50; % Minimum bubble trajectory length (um)
% DG = digraph(s,t,r_inverse); % Directed graph generation
velocity_multiplicator = 1; % Not in use for now
v_propagation = NaN;
v_propagation_manual = 25000; % (um/s) Velocity of the pulse. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330793/
std_hingot_velocity = 0;
debug_propagation_factor = 1; % Propagation slowdown factor
n_paquets = n_bubbles/bb_per_paquet;
%%% Waitbar
f = waitbar(0,'1','Name','Simulation...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);
%%%
for pqt = 1:n_paquets % each paquet of trajectories
    clear bubbles_pqt % bubbles paquet
    bubbles_pqt = cell(bb_per_paquet,1);
    % Check for clicked Cancel button
    for trj = 1:bb_per_paquet % each trajectory
        tic
        if getappdata(f,'canceling')
            delete(f);
            error('Simulation cancelled.')
        end
        %bubbles{jj}.poiseuille = 2;
        %while(bubbles{jj}.poiseuille>1);bubbles{jj}.poiseuille = abs(std*randn(1));end % normal distribution of mean 0 and std = 0.5
        bubbles_pqt{trj}.poiseuille_original = v_poiseuille(floor(length(v_poiseuille)*rand)+1);
        bubbles_pqt{trj}.min_poiseuille_reached = 0;
        bubbles_pqt{trj}.dt = dt;
        if(bubbles_pqt{trj}.poiseuille_original < min_poiseuille) % if poiseuille ratio is lower than threshold
            bubbles_pqt{trj}.poiseuille = min_poiseuille;
            bubbles_pqt{trj}.min_poiseuille_reached = 1;
        else
            bubbles_pqt{trj}.poiseuille = bubbles_pqt{trj}.poiseuille_original;
        end
        %inter_distance = v*dt*(1-bubbles_pqt{trj}.poiseuille); %(um) node_distance with compensation for radial position from center
        clear X Y Z points new_distances dd distances_point_previous distances_next_previous closest_nodes delta pp ax bx cx dx ay by cy dy az bz cz dz  
        %fprintf('2\n');
        while 1 % create new trajectory while too short
            if bypass_N_vs_d_stats == 0
                random_end_node = 1+round(pdfrnd(0:numel(end_nodes_biff_sorted)-1, N_traject_norm, 1));%randi([1 length(end_nodes)],1,1);
            else % bypass
                random_end_node = randi([1 length(end_nodes_biff_sorted)],1,1);
            end
            start = 1;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
            [SP, ~] = shortestpathtree(DG,start,end_nodes_biff_sorted(random_end_node)); % Shortest path
            edges = table2array(SP.Edges);
            nodes = [edges(:,1);edges(end,2)]; % It is the previous node!
            trajectory = pos(nodes,:); % Nodes positions attribution
            d_trajectory = sum(sqrt(sum(diff(trajectory,[],1).^2,2))); % Total length
            x = rand(1); % random distribution
            if and((d_trajectory > min_length),x<bubbles_pqt{trj}.poiseuille)
                break;
            end
        end

        vf_array = [];
        ecg_array = [];
        bubbles_pqt{trj}.d_trajectory = d_trajectory;
        distances = sqrt(sum(diff(trajectory,[],1).^2,2)); % Distances between each original node of the generated trajectory
        distances_cum = cumsum(distances); % Cumulated distances
        xyz = trajectory';
        spline_f = cscvn(xyz); % Creation of the cubic splines
        coefficients = spline_f.coefs; % Getting the coefficients
        start = 0; % (um) % Starting distance from first node
        %%% Vary the distance according to the closest node's radius
        dd = 0;
        new_distances(1,1) = start;
        k = 2;
        bubble_can_go_through = 1;
        while((d_trajectory-new_distances(k-1,1) > max(v_sample_um)*dt)&&bubble_can_go_through)
            previous_nodes_idx = find(distances_cum <= dd); % Finding the nodes that are before the point to get the closest but before node
            if(~isempty(previous_nodes_idx)) % if the starting distance is greater than the first node
                previous_node_idxes(k-1,1) = previous_nodes_idx(end)+1; % This is the previous node's index. 
                distances_point_previous(k-1,1) = dd - distances_cum(previous_node_idxes(k-1,1)-1); % Calculating the distance from previous node to know which of the next and previous are closest
                distances_next_previous(k-1,1) = distances_cum(previous_node_idxes(k-1,1)) - distances_cum(previous_node_idxes(k-1,1)-1);
                if((distances_point_previous(k-1)/distances_next_previous(k-1))<=0.5) % if previous point is closest
                    if(dd<d_trajectory) % if length not exceeding path length
                        closest_nodes(k-1,1) = previous_node_idxes(k-1,1);
                        if(pulsatility==1)
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            wave_delays(k-1) = wave_delay;
                            wave_index = k+period/dt-wave_delay;%+(period/dt)-floor((dd/d_trajectory)*(period/dt));
                            wave_indexes(k-1) = wave_index;
                            vf = v*ecg_normalized(wave_index);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = ecg_normalized(wave_index); % Save ecg_normalized
                            dd = dd + dt*vf;
                        else
                            vf = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = 1; % Save ecg_normalized
                            dd = dd+dt*vf;
                        end
                        new_distances(k,1) = dd;%v_sample(closest_nodes(k-1));%dd + inter_distance*r_norm(closest_nodes(k-1)); % This is the important array which contains the distances between the new nodes
                        k = k+1;
                        if(r(nodes(closest_nodes(k-2))) - padding_bubble)<=0
                            bubble_can_go_through=0;
                        end
                    end
                else % if next node is closest
                    if(dd<d_trajectory) % if length not exceeding path length
                        closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
                        if(pulsatility==1)
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            wave_delays(k-1) = wave_delay;
                            wave_index = k+period/dt-wave_delay;% + (period/dt)-floor((dd/d_trajectory)*(period/dt));
                            wave_indexes(k-1) = wave_index;
                            vf = v*ecg_normalized(wave_index);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = ecg_normalized(wave_index); % Save ecg_normalized
                            dd = dd + dt*vf;
                        else
                            vf = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = 1; % Save ecg_normalized
                            dd = dd + dt*vf;
                        end
                        new_distances(k,1) = dd;%dd + inter_distance*r_norm(closest_nodes(k-1));
                        k = k+1;
                        if(r(nodes(closest_nodes(k-2))) - padding_bubble)<=0
                            bubble_can_go_through=0;
                        end
                    end
                end
            else % the starting distance is less than the first node
                previous_node_idxes(k-1,1) = 1;
                next_node_idx = previous_node_idxes(k-1,1) + 1;
                distances_point_previous(k-1,1) = dd;
                distances_next_previous(k-1,1) = distances_cum(previous_node_idxes(k-1,1));
                if((distances_point_previous(k-1)/distances_next_previous(k-1))<=0.5) % if previous point is closest
                    if(dd<d_trajectory) % if length not exceeding path length
                        closest_nodes(k-1,1) = previous_node_idxes(k-1,1);
                        if(pulsatility==1)
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            wave_index = k+period/dt-wave_delay;% + (period/dt)-floor((dd/d_trajectory)*(period/dt));
                            wave_delays(k-1) = wave_delay;
                            wave_indexes(k-1) = wave_index;
                            vf = v*ecg_normalized(wave_index);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = ecg_normalized(wave_index); % Save ecg_normalized
                            dd = dd + dt*vf;
                        else
                            
                            vf = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = 1; % Save ecg_normalized
                            dd = dd + dt*vf;
                        end
                        new_distances(k,1) = dd;%v_sample(closest_nodes(k-1));%dd + inter_distance*r_norm(closest_nodes(k-1));
                        k = k+1;
                    end
                else
                    if(dd<d_trajectory) % if length not exceeding path length
                        closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
                        if(pulsatility==1)
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            wave_index = k+period/dt-wave_delay;%+ (period/dt)-floor((dd/d_trajectory)*(period/dt)); % The propagation wave
                            wave_delays(k-1) = wave_delay;
                            wave_indexes(k-1) = wave_index;
                            vf = v*ecg_normalized(wave_index);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = ecg_normalized(wave_index); % Save ecg_normalized
                            dd = dd + dt*vf;
                        else
                            vf = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles_pqt{trj}.poiseuille);
                            vf_array(1,k-1) = vf; % save velocity
                            ecg_array(1,k-1) = 1; % Save ecg_normalized
                            dd = dd + dt*vf;
                        end
                        new_distances(k,1) = dd;%v_sample(closest_nodes(k-1));%dd + inter_distance*r_norm(closest_nodes(k-1));
                        k = k+1;
                    end
                end
            end
            dd = new_distances(k-1,1);
        end
%             new_distances = new_distances;
        vf_array(1) = []; % removing first velocity (= 0)
        ecg_array(1) = []; % removing first velocity (= 0)
        bubbles_pqt{trj}.vf_array = vf_array; % saving velocity
        bubbles_pqt{trj}.ecg_array = ecg_array; % saving velocity
        bubbles_pqt{trj}.closest_nodes = nodes(closest_nodes); % Saving the closest nodes indexes
        %%%%% Calculation of the new positions using the cubic splines
        %%%%% coefficients
        L = length(new_distances)-1;
        d = new_distances;
        delta = distances_point_previous./sqrt(distances_next_previous); % distance_point_previous_normalized with the square root. Delta is the scalar used to calculate the position of the new nodes using the distance and the cubic spline
        % point calculation using spline
        pp = (previous_node_idxes(1:L)-1)*3 +1; % Array created to get the good indices of oefficients
        ax = coefficients(pp,1);
        bx = coefficients(pp,2);
        cx = coefficients(pp,3);
        dx = coefficients(pp,4);
        X = ax.*(delta.^3) + bx.*(delta.^2) +...
            cx.*(delta) + dx; % X component
        ay = coefficients(pp+1,1);
        by = coefficients(pp+1,2);
        cy = coefficients(pp+1,3);
        dy = coefficients(pp+1,4);
        Y = ay.*(delta.^3) + by.*(delta.^2) +...
            cy.*(delta) + dy; % Y component
        az = coefficients(pp+2,1);
        bz = coefficients(pp+2,2);
        cz = coefficients(pp+2,3);
        dz = coefficients(pp+2,4);
        Z = az.*(delta.^3) + bz.*(delta.^2) +...
            cz.*(delta) + dz;  % Z component
        XYZ_centerLine = horzcat(X,Y,Z);
        %%% Laminar flow calculation
        clear xyz parallel perpendicular perpendicular2 radii
        parallel = [XYZ_centerLine(2:end,1)-XYZ_centerLine(1:end-1,1) ...
                    XYZ_centerLine(2:end,2)-XYZ_centerLine(1:end-1,2) ...
                    XYZ_centerLine(2:end,3)-XYZ_centerLine(1:end-1,3)]; % vectors parallel to the nodes
        parallel_smooth = smooth(parallel,0.02);
        parallel_smooth = reshape(parallel_smooth,[size(parallel,1) 3]);
        parallel = parallel_smooth;
        perpendicular = zeros(size(parallel,1),3);
        perpendicular2 = zeros(size(parallel,1),3);
        for i = 1:size(parallel,1) % Iterate because vectorised function that gives orthogonal vectors of arrays not found
            perpendiculars = null(parallel(i,:)); % The null() function returns 2 orthogonal vectors to the set of 2 points
            perpendicular(i,:) = perpendiculars(:,1)'; %perpendicular vector 1
            perpendicular2(i,:) = perpendiculars(:,2)'; %perpendicular vector 2
        end
        %%% linear combination
        random_combination1 = rand(1);
        random_combination2 = rand(1);
        lin_combination = (-1+2*random_combination1)*perpendicular+...
                          (-1+2*random_combination2)*perpendicular2;
        %%% Normalize the lin_combination vector to obtain a circular
        %%% distribution rather than a rectangular one
        lin_combination = lin_combination./norm(max(lin_combination));
        %%% Compensate for Poiseuille
        lin_combination = lin_combination.*(sqrt((1-bubbles_pqt{trj}.poiseuille_original))); % compensation of the radial component(lin_combination) by the poiseuille value
    %     lin_combination_smooth = smooth(lin_combination,'loess');
    %     lin_combination_smooth = reshape(lin_combination_smooth,[size(lin_combination,1) 3]);
    %     lin_combination = lin_combination_smooth;
        bubbles_pqt{trj}.radii = abs(r(nodes(closest_nodes)) - padding_bubble); % Radii of the new nodes with compensation with half the bubble size
        %compensation_radii = radii(2:end,1)-radii(1:end-1,1); % Difference in radius between each node
        %compensation_radii = (radii(3:end,1)-radii(1:end-2,1))/2; % Centered difference in radius between each node with error of order 2
        %compensation_radii = smooth(compensation_radii,0.1); % smoothing the derivative curve
        %compensation_rayon_cumul = cumsum(compensation_radii); % Cumulative difference of the radii between each node
        %laminar_xyz = zeros(length(xyz),3);
        %laminar_xyz(1,:) = xyz(1,:)+ lin_combination(1,:).*radii(1); % Point_coordinates + perpendicular_component * radius
        laminar_xyz = XYZ_centerLine(1:end-1,:) + lin_combination.*bubbles_pqt{trj}.radii(1:end-1);
        bubbles_pqt{trj}.XYZ_laminar = laminar_xyz;
        bubbles_pqt{trj}.ID = (pqt-1)*bb_per_paquet + trj;
        %vertices{jj} = laminar_xyz;
        [tot_toc, estimated_time_hours] = DisplayEstimatedTimeOfLoop(tot_toc+toc, bubbles_pqt{trj}.ID, n_bubbles); % Show progress to the user
        time_est_str = ['Estimated time to finish (HH:MM:SS): ' ...
            datestr(estimated_time_hours, 'HH:MM:SS') ' ' ...
            num2str(round(bubbles_pqt{trj}.ID*100/n_bubbles)) '% - ', ...
            num2str(bubbles_pqt{trj}.ID)];
        waitbar(bubbles_pqt{trj}.ID/n_bubbles,f,time_est_str);
    end
    %% Save paquet of bubbles
    save([save_path 'temp\' 'bubbles_pqt_',file_name,'_paquet_',num2str(pqt),'_.mat'],'bubbles_pqt','-v7.3');
    clear bubbles_pqt
end
delete(f)
beep2
%% Gather all Microbubbles
disp('Gathering paquets...')
clear bubbles
bubbles = cell(n_bubbles,1); % Initialization
for pqt = 1:n_paquets
    load([save_path 'temp\' 'bubbles_pqt_',file_name,'_paquet_',num2str(pqt),'_.mat'])
    for trj = 1:bb_per_paquet
        bubbles{trj+(pqt-1)*bb_per_paquet} = bubbles_pqt{trj};
    end
end
clear bubbles_pqt
delete([save_path 'temp\' 'bubbles_pqt_',file_name '*'])
%% Plot trajectories
n_bubbles_plot = 200;
if or(or(display==1,display==2),display==4)
    disp('Plotting Trajectories...');
    figure(6)
    clf
    scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled') % Shortest path nodes);
    n_bubbles = size(bubbles,1);
    for jj = 1:n_bubbles_plot
        hold on
        plot1 = plot3(bubbles{jj}.XYZ_laminar(:,1),bubbles{jj}.XYZ_laminar(:,2),bubbles{jj}.XYZ_laminar(:,3),'LineWidth',1,'Color', [(bubbles{jj}.poiseuille), 0, 1-bubbles{jj}.poiseuille]);
        plot1.Color(4) = 0.3;
%         drawnow
    end
    titre1 = 'Laminar flow simulation with';
    titre2 = num2str(n_bubbles);
    titre3 = ' trajectories.';
    titre_final = [titre1 titre2 ' ' titre3];
    title(titre_final);
    legend('Original nodes','Generated Laminar Flow - Red(fast) Blue(slow)','location','best');
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    zlabel('z','FontSize',20);
    axis equal tight
    view(-22,-22)
end
drawnow
%% Sorting Trajectories as a Function of Flow/Radius
disp('Sorting trajectories...');
clear flow_array
flow_array = [];
stats.RADII = [];
radii_idx = 1;
for ii = 1:n_bubbles % Sorting as a function of r
    flow_array(ii,1) = mean(bubbles{ii}.radii);%bubbles{ii}.d_trajectory*mean(bubbles{ii}.radii)^2 ...
                     %* bubbles{ii}.poiseuille;
    stats.RADII(radii_idx:radii_idx+numel(bubbles{ii}.radii)-1) = bubbles{ii}.radii;
    radii_idx = radii_idx + numel(bubbles{ii}.radii);
end
% flow_array = flow_array/max(flow_array);
[flow_array_sorted,flow_array_idx] = sort(flow_array,1,'descend');
% bubbles_tmp = cell(n_bubbles,1);
% for ii = 1:n_bubbles
%     bubbles_tmp{ii} = bubbles{flow_array_idx(ii)};
% end
% bubbles = bubbles_tmp;
% clear bubbles_tmp
r_mean_sample = linspace(flow_array_sorted(1),flow_array_sorted(end),n_bubbles);
d_mean_sample_log = log(2*r_mean_sample);
N_mean_sample_log = 3.7*d_mean_sample_log -8.2;
N_mean_sample = exp(N_mean_sample_log);
rand_pdf = floor(randpdf(N_mean_sample,1:n_bubbles,[n_bubbles,1]))+1;
rand_pdf_times_N = rand_pdf.*r_mean_sample';
rand_pdf_times_N = floor(rand_pdf_times_N./max(rand_pdf_times_N)...
                    .*max(rand_pdf))+1; % Contains indexes of the bubbles 
rand_pdf_times_N(rand_pdf_times_N > n_bubbles) = n_bubbles; % Fix bound = n_bubbles
%%                                      % to take in the SS calculation
%% Stats
disp('Computing statistics...');
stats.max_d = ceil(max(stats.RADII*2));
stats.min_d = floor(min(stats.RADII*2));
[stats.N_hist,stats.DIAMETER_hist] = hist(stats.RADII*2,(stats.max_d-stats.min_d)/2);
stats.not_zeros_in_N = not(stats.N_hist==0);
stats.N_hist = stats.N_hist(stats.not_zeros_in_N);
stats.DIAMETER_hist = stats.DIAMETER_hist(stats.not_zeros_in_N);
if display == 1
    figure(94);clf;
    hist(stats.RADII*2,10);title('DIAMETERS');
    xlabel('d (\mum)');ylabel('N');
end
stats.x = log(stats.DIAMETER_hist');
stats.y = log(stats.N_hist');
stats.X = [ones(length(stats.x),1) stats.x];
stats.b = stats.X\stats.y;
stats.yCalc2 = stats.X*stats.b;
stats.yHingot = 3.7*stats.x-8.2;
if display == 1
    figure(95);clf;
    scatter(log(stats.DIAMETER_hist),log(stats.N_hist));hold on;plot(stats.x,stats.yCalc2,'--');
    plot(stats.x,stats.yHingot,'*-')
    xlabel('Log diameter');ylabel('Log N');title('Log-Log N vs d')
    legend_title = ['y = ' num2str(stats.b(2)) 'x + ' num2str(stats.b(1))];
    legend('Count per diameter',legend_title,'Ref: y = 3.5x -8.5','Location','Best');
end
%%
if display == 3
    figure(96);clf
    hist(rand_pdf_times_N,100);title('SS Flow Bubbles Probability');
    xlabel('Bubble ID');ylabel('N');
    figure(97);clf
    hist(rand_pdf,100);title('SS Flow Bubble Probability function');
    xlabel('Bubble ID');ylabel('N');
    figure(98);clf
    plot(r_mean_sample,N_mean_sample);
    xlabel('r  mean sample');ylabel('N mean sample');
    figure(99);clf
    subplot(1,2,1);
    n = ceil(max(r_mean_sample));
    hist(flow_array_sorted,n); xlabel('r');ylabel('');
    title('Hist flow array sorted');
    subplot(1,2,2);
    plot(flow_array_sorted)
end
%% Save
disp('Saving bubbles...')
save_file_name = [file_name, '_', num2str(n_bubbles), '_bubbles_', ...
    num2str(n_bubbles_steady_state),'_bubbles_per_frame_', num2str(samp_freq),...
    '_Hz_', num2str(t_steady_state*1000), '_ms_'];
save([save_path 'bubbles_',save_file_name,'.mat'],'bubbles','samp_freq',...
    'n_bubbles','n_bubbles_steady_state','t_steady_state','bubble_size',...
    'pulsatility','slope','intercept','stats','v_propagation_manual',...
    'filename','-v7.3');
%% Steady state flow calculation
disp('Starting steady state flow calculation...');
clear frames frames_velocities
dt = bubbles{1}.dt;
n_frames = t_steady_state/dt;
max_frames = 0;
bubble_count = 0;
loop_counter = 0;
k = 1;  % Bubble per frame counter (row in frames matrix)
ii = 1; % Frame number times 5
frm = 1; % frame number
tot_toc = 0;
% frames_velocities(:,1) = zeros(n_bubbles_steady_state,1);
frames_velocities = NaN(n_bubbles_steady_state,n_frames); % bubbles velocities
frames_ecg = NaN(n_bubbles_steady_state,n_frames); % bubbles ecg amplitude
frames_radii = zeros(n_bubbles_steady_state,n_frames);
frames_poiseuille = zeros(n_bubbles_steady_state,n_frames);
probability_fnct = (linspace(0,1,n_bubbles).^2);
% Simulation time verification
for jj = 1:length(bubbles)
    if size(bubbles{jj}.XYZ_laminar,1) > max_frames
        max_frames = size(bubbles{jj}.XYZ_laminar,1);
    end
end
if 0%n_frames > max_frames
    disp(['The maximum simulation time given the data is : ', num2str(max_frames*bubbles{1}.dt), ' s']);
else
    frames = NaN(n_bubbles_steady_state,n_frames);
    % Generate first set of bubbles
    pp = 1;
    while pp <= n_bubbles_steady_state
        frames(pp,ii) = rand_pdf_times_N(pp); % IDs
        random_index = randi([1 size(bubbles{frames(pp,ii)}.XYZ_laminar,1)],1,1);
        frames(pp,ii+1) = round(floor(random_index/(period/dt))*period/dt +1);    % idx
        if(size(bubbles{frames(pp,ii)}.XYZ_laminar,1)>=frames(pp,ii+1)) 
            frames(pp,ii+(2:4)) = bubbles{frames(pp,ii)}.XYZ_laminar(frames(pp,ii+1),:);
            frames_velocities(pp,1) = bubbles{frames(pp,ii)}.vf_array(frames(pp,ii+1));
            frames_ecg(pp,1) = bubbles{frames(pp,ii)}.ecg_array(frames(pp,ii+1));
            frames_radii(pp,1) = bubbles{frames(pp,ii)}.radii(frames(pp,ii+1));
            frames_poiseuille(pp,1) = bubbles{frames(pp,ii)}.poiseuille;
            pp = pp + 1;
        end
    end
    bubble_count = bubble_count + n_bubbles_steady_state;
    ii = ii + 5;
    frm = frm + 1;
    while frm <= n_frames
        tic
        loop_counter = loop_counter+1;
        while k <= n_bubbles_steady_state % Fill the column with bubbles IDs and time index ii (3 equivalent to 3*dt)
            if(size(bubbles{frames(k,ii-5)}.XYZ_laminar,1) > frames(k,ii-4)) % The trajectory is not ended
                frames(k,ii) = frames(k,ii-5);
                frames(k,ii+1) = frames(k,ii-4)+1;
                frames(k,ii+(2:4)) = bubbles{frames(k,ii)}.XYZ_laminar(frames(k,ii+1),:);
                frames_velocities(k,frm) = bubbles{frames(k,ii)}.vf_array(frames(k,ii+1));
                frames_ecg(k,frm) = bubbles{frames(k,ii)}.ecg_array(frames(k,ii+1));
                frames_radii(k,frm) = bubbles{frames(k,ii)}.radii(frames(k,ii+1));
                frames_poiseuille(k,frm) = bubbles{frames(k,ii)}.poiseuille;
%                 xyz = [frames(k,(ii-5)+(2:4)); frames(k,[ii+(2:4)])];
%                 frames_velocities(k,loop_counter+1) = sqrt(sum(diff(xyz,[],1).^2,2));
%                 frames_velocities(k,loop_counter+1) = smooth(frames_velocities(k,loop_counter+1));
            else
                bubble_count = bubble_count + 1; % add new bubble
                frames(k,ii) = rand_pdf_times_N(bubble_count);
                if(pulsatility == 1)
%                     if(k <= size(bubbles{bubble_count}.wave_delays,2))
%                         sync_pos = period/dt + mod(bubbles{bubble_count}.wave_delays(k),period/dt)+1;%mod(k,period/dt)+1;
%                     else
%                         sync_pos = period/dt + mod(k,period/dt)+1;%mod(k,period/dt)+1;
%                     end
%                     if sync_pos <= size(bubbles{bubble_count}.new_distances,1)
%                         dd = bubbles{bubble_count}.new_distances(sync_pos);
%                     else
%                         dd = 1;
%                     end
%                     wave_index = period/dt-mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                    sync_pos = mod(loop_counter,period/dt)+1;
                    if(sync_pos <= size(bubbles{frames(k,ii)}.XYZ_laminar,1)) % if synchronized position is possible
                        frames(k,ii+1) = sync_pos; % generate position synchronized with frame
                    else
                        frames(k,ii+1) = size(bubbles{frames(k,ii)}.XYZ_laminar,1); % generate position of new bubble at last position
                    end
                else
                    frames(k,ii+1) = randi([1 size(bubbles{frames(k,ii)}.XYZ_laminar,1)],1,1); % generate random position of new bubble
                end
                frames(k,ii+(2:4)) = bubbles{frames(k,ii)}.XYZ_laminar(frames(k,ii+1),:);
                frames_velocities(k,frm) = bubbles{frames(k,ii)}.vf_array(frames(k,ii+1));
                frames_ecg(k,frm) = bubbles{frames(k,ii)}.ecg_array(frames(k,ii+1));
                frames_radii(k,frm) = bubbles{frames(k,ii)}.radii(frames(k,ii+1));
                frames_poiseuille(k,frm) = bubbles{frames(k,ii)}.poiseuille;
            end
            k = k +1;
        end
        ii = ii + 5;
        frm = frm + 1;
        k = 1;
    tot_toc = DisplayEstimatedTimeOfLoop(tot_toc+toc, loop_counter, n_frames);
    end
%     max_velocity = max(max(frames_velocities));
%     frames_velocities = frames_velocities./max_velocity;
    frames_label = ['Bubble ID | Bubble index | X(um) | Y(um) | Z(um)'];
    frames_param.dt = dt;
    frames_param.pulsatility = pulsatility;
    frames_param.t_f = t_f;
    frames_param.n_frames = n_frames;
    save([save_path 'frames_',save_file_name,'.mat'],'frames_label','frames',...
        'frames_velocities','samp_freq','n_bubbles','n_bubbles_steady_state',...
        't_steady_state','bubble_size','pulsatility','filename',...
        'stats','frames_radii','frames_poiseuille','frames_ecg','-v7.3');
end
beep2
disp('Successfully saved frames!')
fprintf('Theoretically, you could have simulated a maximum of  %3.1f s\n',...
    (n_bubbles/bubble_count*t_steady_state));
%% Plot steady state flow
if or(display==1,display==2)
    n_bubbles_steady_state = size(frames,1);
    figure(5);clf
    grid on
    n_frames = size(frames,2)/5;
    fast_forward = 8;
    moving_average = 10;
    for jj = 11:fast_forward:n_frames-moving_average % Starting at 2 since initial velocities are 0
        view_idx = jj/30;
        pp = 3 + (jj-1)*5;
        c = jet(1001);
        scatter3(frames(:,pp),frames(:,pp+2),frames(:,pp+1),3,...
        c(ceil(1000*sqrt(((frames_velocities(:,jj)./max(max(frames_velocities(:,jj-moving_average:jj+moving_average)))))))+1,:));
        axis equal
        xlim([min(pos(:,1)) max(pos(:,1))]); xlabel('x (\mum)');
        ylim([min(pos(:,3)) max(pos(:,3))]); zlabel('z (\mum)');
        zlim([min(pos(:,2)) max(pos(:,2))]); ylabel('y (\mum)');
%         view(54-view_idx*100/samp_freq,21) % For rotating view
        view(135,155);camorbit(180,180)
        set(gca,'GridAlpha',0.5);   
        title([num2str(round(jj/samp_freq,2)) ' s'])
        darkBackground(gcf)
        drawnow
    end
end
%% Plot scatter with speeds
if display == 2
    % Velocities calculation
    max_d = 0;
    dt = bubbles{1}.dt;
    n_bubbles = 3000;%size(bubbles,1);
    for jj = 1:n_bubbles
        if(~isempty(bubbles{jj}.XYZ_laminar))
            if(size(bubbles{jj}.XYZ_laminar,1)>=2)
                %difference = diff(bubbles{jj}.XYZ_laminar,[],1);
                fx_plus_h = bubbles{jj}.XYZ_laminar(3:end,:);
                fx_minus_h = bubbles{jj}.XYZ_laminar(1:end-2,:);
                difference = (fx_plus_h - fx_minus_h)/2; % Centered numerical differentiation of order 2
                bubbles{jj}.velocities = sqrt(sum(difference.^2,2)); % velocity calculation
                if(max(bubbles{jj}.velocities)>max_d)
                    max_d = max(bubbles{jj}.velocities);
                end
                bubbles{jj}.velocities(end+1) = bubbles{jj}.velocities(end);% add a component to have equal length as XYZ
                bubbles{jj}.velocities(end+1) = bubbles{jj}.velocities(end);% add a component to have equal length as XYZ
                bubbles{jj}.velocities = smooth(bubbles{jj}.velocities);
                bubbles{jj}.velocities_normalized = bubbles{jj}.velocities./max_d;
            end
        end
    end
    figure(31);
    clf
    set(gcf,'color','w');
    colormap jet
    for jj = 1:n_bubbles
        if(~isempty(bubbles{jj}.XYZ_laminar))
            if(size(bubbles{jj}.XYZ_laminar,1)>=2)
                n = length(bubbles{jj}.velocities);
                h = scatter3(bubbles{jj}.XYZ_laminar(:,1),...
                    bubbles{jj}.XYZ_laminar(:,2), ...
                    bubbles{jj}.XYZ_laminar(:,3),1,...
                    [bubbles{jj}.velocities_normalized],'Filled');
                alpha = bubbles{jj}.poiseuille;
                set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
                hold on
            end
        end 
    end
    darkBackground(gcf)
    xlabel('\color{white}x (\mum)','FontSize',14);
    ylabel('y (\mum)','FontSize',14);
    zlabel('z (\mum)','FontSize',14);
    title(['\color{white} Bubbles velocities ' num2str(1/dt) 'Hz']);
    c = colorbar;
    c.Color = [1 1 1];
view(-24,46)
end

%% Plot single trajectories in time
if display == 2
    figure(8);
    clf
    set(gcf,'color','w');
    title('Simulation bubbles following a laminar flow');
    hold on
    scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled') % Shortest path nodes);
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    zlabel('z','FontSize',20);
    view(-105,20)
    n_bubbles = size(bubbles,1);
    for jj = 1:n_bubbles
        if(~isempty(bubbles{jj}))
            if(~isempty(bubbles{jj}.XYZ_laminar));plot2 = plot3(bubbles{jj}.XYZ_laminar(:,1), bubbles{jj}.XYZ_laminar(:,2), bubbles{jj}.XYZ_laminar(:,3),'Color', [(bubbles{jj}.poiseuille), 0, 1-bubbles{jj}.poiseuille]);
                plot2.Color(4) = 0.4;end
            drawnow
        end
    end
    legend('Original nodes','Generated trajectories (Red = Fast) (Blue = Slow)');
end
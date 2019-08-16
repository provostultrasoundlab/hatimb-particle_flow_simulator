clear all
close all
clc
%% Path
addpath(genpath('C:\Users\hatim\OneDrive - polymtl.ca\Documents\LabUsons\CODES_HATIM\Rafat\ForHatim'));
%% Visualize
display = 0;
%% Loading
name = 'tree5';
filename = [name '.swc'];
g = importdata(filename);
%% Variables
target = g(:,1);      % IDs
source = g(:,7);  % Parent node
pos = g(:,3:5);   % Positions [x,y,z]
r = g(:,6);%r_temp = g(:,6);       % Radius
r_norm = r./max(r);
r_inverse = 1./r;
r_inverse_norm = 1./r_norm;
%% Positions scaling
% dim_1 = 2; % (um)
% dim_2 = 2; % (um)
% dim_3 = 50; % (um)
% pos(:,1) = pos(:,1) * dim_1;
% pos(:,2) = pos(:,2) * dim_2;
% pos(:,3) = pos(:,3) * dim_3;
% r = r .* (dim_1+dim_2)/2;
%% Smoothing radius
r_orig_smooth = smooth(r,0.1);
r_smooth = smooth(r_norm,0.1); % Smoothed Radius
r_norm = r_smooth;
%% Affichage rayons
if display == 1
    figure;plot(r);hold on;plot(r_orig_smooth(1:end-10),'LineWidth',5);
    legend('r','r smooth');
    title('Rayon en fonction du noeud du reseau');
end
%% Graph
s = source+2;
t = target+2;
C = setxor(s,1:length(s));
end_nodes = [C-1;s(end)];
clear C
%% Finding bifurcations
clear A
A = s;
[uniqueA i j] = unique(A,'first');
indexToDupes = find(not(ismember(1:numel(A),i)));
biff_nodes = s(indexToDupes-1)-1;
% figure(72);
% clf
% %scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled');
% plot3(pos(:,1),pos(:,2),pos(:,3),'.k');
% hold on
% plot3(pos(biff_nodes,1),pos(biff_nodes,2),pos(biff_nodes,3),'og') % Starting flow node
% title('Bifurcations');
% %% Directed graph
% clf
% start = 2;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
% finish = start + randi([1 size(s,1)-start-1],1,1); % random integer from [start:#target_nodes-start]
% DG = digraph(s,t); % Directed graph generation
% [SP, D] = shortestpathtree(DG,start,finish); % Shortest path
% edges = table2array(SP.Edges); % Conversion
% nodes = [edges(:,1);edges(end,2)];
% trajectory = pos(nodes-1,:);
% %Plot 1 trajectory
% figure(1);f = plot(SP); % Plot the shortest path graph variable
% f.XData = [pos(1,1);pos(:,1)];
% f.YData = [pos(1,2);pos(:,2)];
% f.ZData = [pos(1,3);pos(:,3)];
% f.EdgeColor = 'r';
% f.LineWidth = 3;
% f.ArrowSize = 10;
% hold on
% plot3(pos(biff_nodes,1),pos(biff_nodes,2),pos(biff_nodes,3),'og') % Starting flow node
% xlabel('x','FontSize',20);
% ylabel('y','FontSize',20);
% zlabel('z','FontSize',20);
% view(20,30)
% set(gcf,'color','w');
%% Volume calculation
total_vessels_length = 0;
vectors_volume_calc = pos(t(1:length(t)-1),:)-pos(s(1:length(t)-1),:); % Parallel vectors between all nodes
distances_volume_calc = sqrt(sum(diff(vectors_volume_calc,[],1).^2,2)); % Euclidian norms calculation
length_volume_calc = sum(distances_volume_calc);% um
mean_radius = mean(r_orig_smooth); % um
total_vessel_volume = length_volume_calc*pi()*mean_radius.^2; % um^3
disp(['Total vasculature volume = ' num2str(total_vessel_volume*1E-9) ' mm^3']);
whole_volume = (max(pos(:,1))-min(pos(:,1)))*(max(pos(:,2))-min(pos(:,2)))*(max(pos(:,3))-min(pos(:,3)));
disp(['Vessel volume ratio = ' num2str(total_vessel_volume/whole_volume*100) '%']);
vessel_volume_mL = total_vessel_volume*1E-9/1000;
C_MB = 2E5; % Concentration of MicroBubbles in the bloodstream, Hingot et al, 2019
MB = round(vessel_volume_mL*C_MB); % Number of microbubbles in the studied vessels
disp(['Number of MB in vessels : ' num2str(MB)]);
%% Affichage des points extremes
if display == 1
    DG_display = digraph(s,t,r_inverse); % Directed graph generation
    figure(99);clf;f = plot(DG_display,'Layout','force');
    f.XData = [pos(1,1);pos(:,1)];
    f.YData = [pos(1,2);pos(:,2)];
    f.ZData = [pos(1,3);pos(:,3)];
    f.EdgeColor = 'k';
    f.LineWidth = 1.5;
    f.ArrowSize = 5;
    view(20,30);
    hold on
    plot3(pos(end_nodes,1),pos(end_nodes,2),pos(end_nodes,3),'og');
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    zlabel('z','FontSize',20);
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
d_sample_log = log(2*r_orig_smooth);
v_sample_log = 1.9*(d_sample_log) -6;
v_sample = exp(v_sample_log);
v_sample_um = v_sample*1000;
d_sample = 2*r_orig_smooth;
N_sample_log = 3.7*d_sample_log -8.2;
N_sample = exp(N_sample_log);
%%% Affichage
if display == 1 
    figure;
    subplot(2,2,1);plot(log_d,log_N,'LineWidth',2);
    hold on; plot(d_sample_log,N_sample_log,'.');
    grid on;title('Dependency of the bubble rate with vessel’s diameter');xlabel('log(d)');ylabel('log(N)');
    axis([0 5 0 8]);
    subplot(2,2,2);plot(d,N,'LineWidth',2);
    hold on; plot(d_sample,N_sample,'.');
    grid on;title('Dependency of the bubble rate with vessel’s diameter');xlabel('d (um)');ylabel('N');
    %%% v
    subplot(2,2,3);plot(log_d,log_v,'LineWidth',2);
    hold on; plot(d_sample_log,v_sample_log,'.');
    legend('Hingot','Sample data','Location','Best');
    grid on;title('Dependency of maximum velocity with vessel’s diameter');xlabel('log(d)');ylabel('log(v)');
    axis([0 5 -5 4]);
    subplot(2,2,4);plot(d,v,'LineWidth',2);
    hold on; plot(d_sample,v_sample,'.');
    legend('Hingot','Sample data','Location','Best');
    grid on;title('Dependency of maximum velocity with vessel’s diameter');xlabel('d (um)');ylabel('v (mm/s)');
    %figure;plot(d,v,'LineWidth',2);xlabel('d (um)');ylabel('v');title('Sample velocities');
    %axis([15 45 0 3.5]);
end
%% Poiseuille
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
samp_freq = 10000; %(Hz)
dt = 1/samp_freq;
t_f = 100;
x = 0:dt:t_f-dt;
ecg_raw = ecg(BPM,dt,t_f);
ecg_filtered = ecg_raw-min(ecg_raw);
ecg_filtered2 = ecg_filtered./max(ecg_filtered);
ecg_filtered3 = ecg_filtered2+0.5;
ecg_normalized = ecg_filtered3;
% figure;plot(x,ecg_normalized);xlabel('time(s)');ylabel('Normalized Amplitude');
% figure(44)
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
%% Plotting ecg in time
if display == 1
    Make_video = 0;
    figure(22)
    clf
    if Make_video==1
        v = VideoWriter(['pulse_',num2str(t_f),'s.mp4']);
        v.Quality = 100;
        v.FrameRate = samp_freq/5;
        open(v);
    end
    for i = 1:length(x)
        plot(x(1:i),ecg_normalized(1:i),'g','LineWidth',2); hold on
        xlabel('Time (s)')
        ylabel('Pulse')
        darkBackground(gcf)
        axis([0+(i-41)*dt 1+(i-41)*dt min(ecg_normalized) max(ecg_normalized)])
        plot(x(i),ecg_normalized(i),'.r','MarkerSize',20)
        set(gca,'FontSize',12)
        pause(dt)
        if Make_video==1
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        clf
    end
    if Make_video==1 ; close(v); end
end
%% Simuation
clear bubbles
n_bubbles = 1000; % Number of bubbles trajectories generated
bubble_size = 5; % Bubble diameter (um)
padding_bubble = bubble_size/2; % To account for the fact that the bubbles are not infinitesimal points
bubbles = cell(n_bubbles,1); % Initialization
tot_toc = 0; % For displaying progress to the user
min_length = 200; % Minimum bubble trajectory length (um)
min_poiseuille = 0.2; % Minimum Poiseuille value (a value of 0 causes an infinite computation time since the bubble doesn't move)
DG = digraph(s,t,r_inverse); % Directed graph generation
pulsatility = 1; % 1 = Yes | 0 = No
velocity_multiplicator = 1; % Multiplies velocities according to Hingot V et al, 2019
v_propagation = NaN;
v_propagation_manual = 88000;
std_hingot_velocity = 0;
debug_propagation_factor = 1; % Propagation slowdown factor
for jj = 1:n_bubbles
    tic
    %bubbles{jj}.poiseuille = 2;
    %while(bubbles{jj}.poiseuille>1);bubbles{jj}.poiseuille = abs(std*randn(1));end % normal distribution of mean 0 and std = 0.5
    bubbles{jj}.poiseuille_original = v_poiseuille(floor(length(v_poiseuille)*rand)+1);
    bubbles{jj}.min_poiseuille_reached = 0;
    bubbles{jj}.dt = dt;
    if(bubbles{jj}.poiseuille_original < min_poiseuille) % if poiseuille ratio is lower than threshold
        bubbles{jj}.poiseuille = min_poiseuille;
        bubbles{jj}.min_poiseuille_reached = 1;
    else
        bubbles{jj}.poiseuille = bubbles{jj}.poiseuille_original;
    end
    %inter_distance = v*dt*(1-bubbles{jj}.poiseuille); %(um) node_distance with compensation for radial position from center
    clear X Y Z points new_distances dd distances_point_previous distances_next_previous closest_nodes delta pp ax bx cx dx ay by cy dy az bz cz dz  
    %fprintf('2\n');
    while 1 % create new trajectory while too short
        random_end_node = randi([1 length(end_nodes)],1,1);
        start = 2;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
        [SP, ~] = shortestpathtree(DG,start,end_nodes(random_end_node)); % Shortest path
        edges = table2array(SP.Edges);
        nodes = [edges(:,1);edges(end,2)]-1; % It is the previous node!
        trajectory = pos(nodes,:); % Nodes positions attribution
        d_trajectory = sum(sqrt(sum(diff(trajectory,[],1).^2,2))); % Total length
        if(d_trajectory > min_length)
            break;
        end
    end
    if(1)
        bubbles{jj}.d_trajectory = d_trajectory;
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
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            bubbles{jj}.wave_delays(k-1) = wave_delay;
                            wave_index = k+period/dt-wave_delay;%+(period/dt)-floor((dd/d_trajectory)*(period/dt));
                            bubbles{jj}.wave_indexes(k-1) = wave_index;
                            dd = dd+dt*v*ecg_normalized(wave_index);
                        else
                            dd = dd+dt*v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                        end
                        new_distances(k,1) = dd;%v_sample(closest_nodes(k-1));%dd + inter_distance*r_norm(closest_nodes(k-1)); % This is the important array which contains the distances between the new nodes
                        k = k+1;
                        if(r_orig_smooth(nodes(closest_nodes(k-2))) - padding_bubble)<=0
                            bubble_can_go_through=0;
                        end
                    end
                else % if next node is closest
                    if(dd<d_trajectory) % if length not exceeding path length
                        closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
                        if(pulsatility==1)
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            bubbles{jj}.wave_delays(k-1) = wave_delay;
                            wave_index = k+period/dt-wave_delay;% + (period/dt)-floor((dd/d_trajectory)*(period/dt));
                            bubbles{jj}.wave_indexes(k-1) = wave_index;
                            dd = dd+dt*v*ecg_normalized(wave_index);
                        else
                            dd = dd+dt*v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                        end
                        new_distances(k,1) = dd;%dd + inter_distance*r_norm(closest_nodes(k-1));
                        k = k+1;
                        if(r_orig_smooth(nodes(closest_nodes(k-2))) - padding_bubble)<=0
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
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                            if k == 2 % Finding maximum velocity at begining of trajectory and set it as propagation velocity
                                v_propagation = v_propagation_manual;%v/debug_propagation_factor;
                            end
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            wave_index = k+period/dt-wave_delay;% + (period/dt)-floor((dd/d_trajectory)*(period/dt));
                            bubbles{jj}.wave_delays(k-1) = wave_delay;
                            bubbles{jj}.wave_indexes(k-1) = wave_index;
                            dd = dd+dt*v*ecg_normalized(wave_index);
                        else
                            dd = dd+dt*v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
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
                            v = v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                            wave_delay = mod(floor((dd/v_propagation)*(period/dt)),period/dt);
                            wave_index = k+period/dt-wave_delay;%+ (period/dt)-floor((dd/d_trajectory)*(period/dt)); % The propagation wave
                            bubbles{jj}.wave_delays(k-1) = wave_delay;
                            bubbles{jj}.wave_indexes(k-1) = wave_index;
                            dd = dd+dt*v*ecg_normalized(wave_index);
                        else
                            dd = dd+dt*v_sample_um(nodes(closest_nodes(k-1)))*velocity_multiplicator*(bubbles{jj}.poiseuille);
                        end
                        new_distances(k,1) = dd;%v_sample(closest_nodes(k-1));%dd + inter_distance*r_norm(closest_nodes(k-1));
                        k = k+1;
                    end
                end
            end
            dd = new_distances(k-1,1);
        end
        bubbles{jj}.new_distances = new_distances;
        bubbles{jj}.closest_nodes = nodes(closest_nodes); % Saving the closest nodes indexes
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
        X = ax.*(delta.^3) + bx.*(delta.^2) + cx.*(delta) + dx; % X component
        ay = coefficients(pp+1,1);
        by = coefficients(pp+1,2);
        cy = coefficients(pp+1,3);
        dy = coefficients(pp+1,4);
        Y = ay.*(delta.^3) + by.*(delta.^2) + cy.*(delta) + dy; % Y component
        az = coefficients(pp+2,1);
        bz = coefficients(pp+2,2);
        cz = coefficients(pp+2,3);
        dz = coefficients(pp+2,4);
        Z = az.*(delta.^3) + bz.*(delta.^2) + cz.*(delta) + dz;  % Z component
        bubbles{jj}.XYZ_centerLine = horzcat(X,Y,Z);
    end
    %%% Laminar flow calculation
    clear xyz parallel perpendicular perpendicular2 radii
    xyz = bubbles{jj}.XYZ_centerLine;
    parallel = [xyz(2:end,1)-xyz(1:end-1,1) xyz(2:end,2)-xyz(1:end-1,2) xyz(2:end,3)-xyz(1:end-1,3)]; % vectors parallel to the nodes
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
    lin_combination = (-1+2*rand(1))*perpendicular+(-1+2*rand(1))*perpendicular2;
    %%% Normalize the lin_combination vector to obtain a circular
    %%% distribution rather than a rectangular one
    lin_combination = lin_combination./norm(max(lin_combination));
    %%% Compensate for Poiseuille
    lin_combination = lin_combination.*(sqrt((1-bubbles{jj}.poiseuille_original))); % compensation of the radial component(lin_combination) by the poiseuille value
%     lin_combination_smooth = smooth(lin_combination,'loess');
%     lin_combination_smooth = reshape(lin_combination_smooth,[size(lin_combination,1) 3]);
%     lin_combination = lin_combination_smooth;
    radii = abs(r_orig_smooth(nodes(closest_nodes)) - padding_bubble); % Radii of the new nodes with compensation with half the bubble size
    bubbles{jj}.radii = radii;
    %compensation_radii = radii(2:end,1)-radii(1:end-1,1); % Difference in radius between each node
    %compensation_radii = (radii(3:end,1)-radii(1:end-2,1))/2; % Centered difference in radius between each node with error of order 2
    %compensation_radii = smooth(compensation_radii,0.1); % smoothing the derivative curve
    %compensation_rayon_cumul = cumsum(compensation_radii); % Cumulative difference of the radii between each node
    %laminar_xyz = zeros(length(xyz),3);
    %laminar_xyz(1,:) = xyz(1,:)+ lin_combination(1,:).*radii(1); % Point_coordinates + perpendicular_component * radius
    laminar_xyz = xyz(1:end-1,:) + lin_combination.*radii(1:end-1);
    bubbles{jj}.XYZ_laminar = laminar_xyz;
    bubbles{jj}.ID = jj;
    %vertices{jj} = laminar_xyz;
    tot_toc = DisplayEstimatedTimeOfLoop(tot_toc+toc, jj, n_bubbles); % Show progress to the user
end
save('bubbles_10KHz_3000bb_Erwan.mat','bubbles','-v7.3');
beep2
%% Steady state flow calculation
clear frames frames_velocities
MB = 100;%size(bubbles,1)/5;
t_f = 10; % Desired simulation time (s)
dt = bubbles{1}.dt;
n_frames = t_f/dt;
max_frames = 0;
bubble_count = 0;
loop_counter = 0;
k = 1;  % Bubble per frame counter
ii = 1; % Frame number
tot_toc = 0;
frames_velocities(:,1) = zeros(MB,1);
% Simulation time verification
for jj = 1:length(bubbles)
    if size(bubbles{jj}.XYZ_laminar,1) > max_frames
        max_frames = size(bubbles{jj}.XYZ_laminar,1);
    end
end
if n_frames > max_frames
    disp(['The maximum simulation time given the data is : ', num2str(max_frames*bubbles{1}.dt), ' s']);
else
    frames = NaN(MB,n_frames);
    % Generate first set of bubbles
    for pp = 1:MB
        frames(pp,ii) = pp; % IDs
        random_index = randi([1 size(bubbles{frames(pp,ii)}.XYZ_laminar,1)],1,1);
        frames(pp,ii+1) = round(floor(random_index/(period/dt))*period/dt +1);    % idx
        frames(pp,ii+2:5) = bubbles{frames(pp,ii)}.XYZ_laminar(frames(pp,ii+1),:);
    end
    bubble_count = bubble_count + MB;
    ii = ii + 5;
    while ii <= n_frames*5
        tic
        loop_counter = loop_counter+1;
        while k <= MB % Fill the column with bubbles IDs and time index ii (3 equivalent to 3*dt)
            if(size(bubbles{frames(k,ii-5)}.XYZ_laminar,1) > frames(k,ii-4)) % The trajectory is not ended
                frames(k,ii) = frames(k,ii-5);
                frames(k,ii+1) = frames(k,ii-4)+1;
                frames(k,ii+(2:4)) = bubbles{frames(k,ii)}.XYZ_laminar(frames(k,ii+1),:);
                xyz = [frames(k,(ii-5)+(2:4)); frames(k,[ii+(2:4)])];
                frames_velocities(k,loop_counter+1) = sqrt(sum(diff(xyz,[],1).^2,2));
                frames_velocities(k,loop_counter+1) = smooth(frames_velocities(k,loop_counter+1));
            else
                bubble_count = bubble_count + 1; % add new bubble
                frames(k,ii) = bubble_count;
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
                frames_velocities(k,loop_counter) = 0;
            end
            k = k +1;
        end
        ii = ii + 5;
        k = 1;
    tot_toc = DisplayEstimatedTimeOfLoop(tot_toc+toc, loop_counter, n_frames);
    end
    max_velocity = max(max(frames_velocities));
    frames_velocities = frames_velocities./max_velocity;
    frames_label = ['Bubble ID | Bubble index | X(um) | Y(um) | Z(um)'];
    frames_param.dt = dt;
    frames_param.pulsatility = pulsatility;
    frames_param.t_f = t_f;
    frames_param.n_frames = n_frames;
    save('frames_10KHz_400bb_10s_Erwan.mat','frames_label','frames','frames_velocities','-v7.3');
end
beep2
%% Plot steady state flow
Make_video = 0;
MB = size(frames,1);
figr66 = figure(66);
grid on
clf
n_frames = size(frames,2)/5;
fast_forward = 1;
set(gcf,'color','w');
colormap jet
if Make_video==1
    v = VideoWriter('bubble_flow_puls_v4.mp4');
    v.Quality = 100;
    v.FrameRate = samp_freq/5;
    open(v);
end
for jj = 1:fast_forward:n_frames-1
    view_idx = jj/30;
    pp = 3 + (jj-1)*5;
%     h = scatter3(pos(:,1),pos(:,2),pos(:,3),1,[1 1 1]); % Shortest path nodes);
%     alpha = 0.2;
%     set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
%     hold on
    %plot3(frames(:,pp),frames(:,pp+1),frames(:,pp+2),'or');
    c = jet(1001);
    scatter3(frames(:,pp),frames(:,pp+2),frames(:,pp+1),15,c(ceil(1000*sqrt(((frames_velocities(:,jj)))))+1,:));%(sqrt(sqrt(frames_velocities(:,jj)))) ; ones(MB,1)
%     ax.gridAlpha = 0.5;
    %     grid off
    axis equal
%     xlim([-5 max(pos(:,1))])
%     ylim([-5 max(pos(:,3))])
%     zlim([-5 max(pos(:,2))])
    xlim([-5 200])
    ylim([-5 175])
    zlim([-5 200])
%     xlabel('\color{white}x (\mum)','FontSize',14);
%     ylabel('z (\mum)','FontSize',14);
%     zlabel('y (\mum)','FontSize',14);
%     if pulsatility ==1
% %         title(['\color{white}Steady state flow simulation with pulsatility : ' num2str(MB) ' bubbles | Time (s): ' num2str(jj*dt,'%4.3f')],'FontSize',20);
%         title(['\color{white}Time (s): ' num2str(jj*dt,'%4.3f')],'FontSize',20);
%     else
%         title(['\color{white}Steady state flow simulation : ' num2str(MB) ' bubbles | Time (s): ' num2str(jj*bubbles{1}.dt,'%4.3f')]);
%     end
    %title(num2str(jj*bubbles{1}.dt));
    view(-45,15)
%     view(54,21)
%     view(156,22)
%     view(-100,5)
%     view(-139+view_idx,24)
%     view(54-view_idx*1000/samp_freq,21)
%     c_bar = colorbar;
%     c_bar.Color = [1 1 1];
%     c_bar.Label.String = 'Normalized Velocities^2';
    size_axis = 50;
    hold on
    axis_matrix_x = [0 0 0;size_axis 0 0];
    axis_matrix_y = [0 0 0;0 size_axis 0];
    axis_matrix_z = [0 0 0;0 0 size_axis];
    quiver3(axis_matrix_x(1,1),axis_matrix_x(1,2),axis_matrix_x(1,3),axis_matrix_x(2,1),axis_matrix_x(2,2),axis_matrix_x(2,3),'r','MaxHeadSize',.5,'linewidth',5)
    quiver3(axis_matrix_y(1,1),axis_matrix_y(1,2),axis_matrix_y(1,3),axis_matrix_y(2,1),axis_matrix_y(2,2),axis_matrix_y(2,3),'g','MaxHeadSize',.5,'linewidth',5)
    quiver3(axis_matrix_z(1,1),axis_matrix_z(1,2),axis_matrix_z(1,3),axis_matrix_z(2,1),axis_matrix_z(2,2),axis_matrix_z(2,3),'b','MaxHeadSize',.5,'linewidth',5)
%     darkBackground(gcf)   
    set(gca,'XColor', 'none','YColor','none','ZColor','none')
%     pause(dt*fast_forward)
    set(gca,'GridAlpha',0.5);    
    pause
    hold off
    if Make_video==1
        frame = getframe(gcf);
        writeVideo(v,frame);  
    end
end
if Make_video==1 ; close(v); end

%% Plot Pulsatility
figure(72)
clf
for p = 1:100%size(frames,2)/5
    ii = (p-1)*5 + 1; % Column index
    for m = 1:size(frames,1)
        plot(m,bubbles{frames(m,ii)}.new_distances(frames(m,ii+1)),'.');
        hold on
    end
    xlabel('Bubbles')
    ylabel('Trajectory length (um)')
    axis([1 size(frames,1) 0 3000])
    pause(0.5)
    hold off
end 
%% Plot frames Velocities in time
figure(32)
clf
for jj = 1:1000%size(frames_velocities,1) % for each bubble
    jj
    n = size(frames_velocities,2);
    x1 = dt:dt:n*dt;
    plot(x1(1:200),frames_velocities(jj,1:200))
    hold on
    %axis([0 1 0 0.3])
    title(jj)
%     pause(0.1)
end
titre = ['Bubble velocity in time - ' num2str(1/dt) 'Hz'];
title(titre)
xlabel('Time (s)')
ylabel('Normalized Velocities');
grid on
set(gca,'FontSize',14)

%% Plot scatter with speeds
tic
% Velocities calculation
max_d = 0;
dt = bubbles{1}.dt;
n_bubbles = 500;%size(bubbles,1);
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
            h = scatter3(bubbles{jj}.XYZ_laminar(:,1), bubbles{jj}.XYZ_laminar(:,2), bubbles{jj}.XYZ_laminar(:,3),1,[bubbles{jj}.velocities_normalized]);
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
toc
figure(32)
clf
for jj = 1:n_bubbles
    n = size(bubbles{jj}.velocities,1);
    x1 = dt:dt:n*dt;
    plot(x1,bubbles{jj}.velocities)
    hold on
end
titre = ['Bubble velocity in time (pulsatility) ' num2str(1/dt) 'Hz'];
title(titre)
xlabel('Time (s)')
ylabel('Normalized Velocity');
grid on
axis([0 3 0 1]);
%saveas(gcf,'laminar-2019-06-13-10000.png');

%% plot dark
figr = figure
plot3(1:10,1:10,1:10,'.')
%% Turn 3d plot in time
for i = 1:3600
    view(-45+i/10,45)
    pause(0.5)
end
%% Plot scatter with speeds
% Velocities calculation
max_d = 0;
n_bubbles = 1000;%size(bubbles,1);
for jj = 1:n_bubbles
    if(~isempty(bubbles{jj}.XYZ_laminar))
        if(size(bubbles{jj}.XYZ_laminar,1)>=2)
            bubbles{jj}.velocities = sqrt(sum(diff(bubbles{jj}.XYZ_laminar,[],1).^2,2));
            if(max(bubbles{jj}.velocities)>max_d)
                max_d = max(bubbles{jj}.velocities);
            end
            bubbles{jj}.velocities = bubbles{jj}.velocities./max_d;
            bubbles{jj}.velocities(end+1) = bubbles{jj}.velocities(end);
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
            h = scatter3(bubbles{jj}.XYZ_laminar(:,1), bubbles{jj}.XYZ_laminar(:,2), bubbles{jj}.XYZ_laminar(:,3),1,[bubbles{jj}.velocities],'filled');
            alpha = bubbles{jj}.poiseuille;
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
            hold on
        end
    end 
end
%saveas(gcf,'laminar-2019-06-13-10000.png');
%% Test Streamparticles
% figure(97);clf;streamparticles(vertices);
%% Dérivées
figure;
plot(compensation_radii); hold on
plot(compensation_radii_centree);
plot(compensation_radii_smooth);
legend('Originale','Ordre 2 centrée','Ordre smooth');
%% Affichage
figure(13)
clf
%scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled') % Shortest path nodes);
n_bubbles = size(bubbles,1);
for jj = 1:n_bubbles
hold on
    plot1 = plot3(bubbles{jj}.XYZ_laminar(:,1),bubbles{jj}.XYZ_laminar(:,2),bubbles{jj}.XYZ_laminar(:,3),'LineWidth',1,'Color', [(bubbles{jj}.poiseuille), 0, 1-bubbles{jj}.poiseuille]);
    plot1.Color(4) = 0.05;
    jj
end
titre1 = 'Laminar flow simulation with ';
titre2 = num2str(n_bubbles);
titre3 = ' trajectories.';
% titre4 = num2str(tot_toc);
% titre5 = ' s.';
titre_final = [titre1 titre2 ' ' titre3];
title(titre_final);
legend('Original nodes','Generated Laminar Flow - Red(fast) Blue(slow)','location','best');
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
view(-75,20)
%% Plot single trajectories in time
figure(8);
clf
set(gcf,'color','w');
title('Simulation of 1000 bubbles following a laminar flow');
hold on
scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled') % Shortest path nodes);
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
view(-105,20)
n_bubbles = size(bubbles,1);
for jj = 1:n_bubbles
    if(~isempty(bubbles{jj}))
        %if(~isempty(bubbles{jj}.distances));plot1 = plot3(bubbles{jj}.gen_points(:,1), bubbles{jj}.gen_points(:,2), bubbles{jj}.gen_points(:,3),'b');plot1.Color(4) = 0.2;end
        if(~isempty(bubbles{jj}.XYZ_laminar));plot2 = plot3(bubbles{jj}.XYZ_laminar(:,1), bubbles{jj}.XYZ_laminar(:,2), bubbles{jj}.XYZ_laminar(:,3),'Color', [(bubbles{jj}.poiseuille), 0, 1-bubbles{jj}.poiseuille]);
            plot2.Color(4) = 0.4;end
        pause(0.1);
    end
    
end
legend('Original nodes','Generated trajectories (Red = Fast) (Blue = Slow)');
% titre = ['Simulation of ',num2str(jj),' bubbles following a laminar flow'];
% title(titre);
%legend('Before','Smoothing using csaps','location','best');
%% Plot velocity and diameter dependancy
% Velocities calculation
max_d = 0;
n_bubbles = 1000;%size(bubbles,1);
for jj = 1:n_bubbles
    jj
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
            bubbles{jj}.velocities = bubbles{jj}.velocities./dt;
            bubbles{jj}.velocities = smooth(bubbles{jj}.velocities);
            bubbles{jj}.velocities_normalized = bubbles{jj}.velocities./max_d;
        end
    end
end
log_d = linspace(0,5,1000);
log_N = 3.7*log_d -8.2;
d = exp(log_d);
N = exp(log_N);
log_v = 1.9*log_d -6;
v = exp(log_v);
d_sample_log = log(2*r_orig_smooth);
v_sample_log = 1.9*(d_sample_log) -6;
v_sample = exp(v_sample_log);
v_sample_um = v_sample*1000;
d_sample = 2*r_orig_smooth;
N_sample_log = 3.7*d_sample_log -8.2;
N_sample = exp(N_sample_log);
figure(21);
clf
set(gcf,'color','w');
colormap jet
loglog(d,v,'Color',[.6 0 0],'LineWidth',5);
hold on
grid on
axis([5 45 0 4]);
for jj = 1:n_bubbles
    if(~isempty(bubbles{jj}.XYZ_laminar))
        jj
%         bubbles{jj}.distances = sqrt(sum(diff(bubbles{jj}.XYZ_laminar,[],1).^2,2));
%         v_plot = bubbles{jj}.distances./dt./1000;
        v_plot = bubbles{jj}.velocities./1000;
        d_plot = 2*r_orig_smooth(bubbles{jj}.closest_nodes(2:end));
        loglog(d_plot,v_plot,'.','MarkerSize',1,'MarkerEdgeColor',[1-bubbles{jj}.poiseuille 1-bubbles{jj}.poiseuille 1-bubbles{jj}.poiseuille])
        x = log(d_plot);
        y = log(v_plot);
        p = polyfit(x,y, 1);
        m(jj) = p(1);
        k(jj) = p(2);
        hold on
    end 
end
xlabel('Diameter (um)');
ylabel('Velocity (mm/s)');
legend('Hingot V, et al.','Gen data', 'location','best');
mean_m = mean(m);
std_m = std(m,1,'all');
mean_k = mean(k);
fprintf('Mean p1:%3.2f, mean p2:%3.2f\n',mean_m,mean_k);
log_v_sample_mean = mean_m*log_d +mean_k;
v_mean = exp(log_v_sample_mean);
% title(['Velocity and Vessel Diameter Dependancy | ' 'Mean: y = ' num2str(mean_m) 'x + ' num2str(mean_k)])
title(['Velocity and Vessel Diameter Dependancy'])
set(gca,'fontweight','bold','FontSize',28)
set(gca,'TickDir','out'); % The only other option is 'in'

%loglog(d,v_mean,'Color',[0.8 0 0],'LineWidth',5); % mean of data

%% Plot velocity and diameter dependancy NEWWW
% Velocities calculation
max_d = 0;
dt = bubbles{1}.dt;
n_bubbles = 100;size(bubbles,1);
for jj = 1:n_bubbles
    jj
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
            bubbles{jj}.velocities = bubbles{jj}.velocities./dt;
            bubbles{jj}.velocities = smooth(bubbles{jj}.velocities);
            bubbles{jj}.velocities_normalized = bubbles{jj}.velocities./max_d;
        end
    end
end
log_d = linspace(0,5,1000);
log_N = 3.7*log_d -8.2;
d = exp(log_d);
N = exp(log_N);
log_v = 1.9*log_d -6;
v = exp(log_v);
d_sample_log = log(2*r_orig_smooth);
v_sample_log = 1.9*(d_sample_log) -6;
v_sample = exp(v_sample_log);
v_sample_um = v_sample*1000;
d_sample = 2*r_orig_smooth;
N_sample_log = 3.7*d_sample_log -8.2;
N_sample = exp(N_sample_log);
figure(23);
clf
set(gcf,'color','w');
colormap jet
loglog(d,v,'Color',[0 0 0],'LineWidth',5);
hold on
grid on
axis([5 45 0 4]);
X = []; Y = [];
for jj = 1:n_bubbles
    if(~isempty(bubbles{jj}.XYZ_laminar))
        jj
%         bubbles{jj}.distances = sqrt(sum(diff(bubbles{jj}.XYZ_laminar,[],1).^2,2));
%         v_plot = bubbles{jj}.distances./dt./1000;
        v_plot = bubbles{jj}.velocities./1000;
        d_plot = 2*r_orig_smooth(bubbles{jj}.closest_nodes(2:end));
        loglog(d_plot,v_plot,'.')
        x = log(d_plot);
        X = vertcat(X,x);
        y = log(v_plot);
        Y = vertcat(Y,y);
        p = polyfit(x,y, 1);
        m(jj) = p(1);
        k(jj) = p(2);
        hold on
    end 
end

%%% Polyfit
P = polyfit(X,Y,1);
P(1)
P(2)
Yfit = polyval(P,X);
Yresid = Y - Yfit;
SSresid = sum(Yresid.^2);
SStotal = (length(Y)-1) * var(Y);
rsq = 1 - SSresid/SStotal;


xlabel('Diameter (um)');
ylabel('Velocity (mm/s)');
axis([5 45 0.01 10]);
yticklabels([0.01 0.1 1 10]);
set(gca,'FontSize',20,'FontWeight','bold','MinorGridColor',[0 0 0])
legend('Hingot V, et al.','Gen data', 'location','best');
fprintf(['Fit: y = ' num2str(P(1)) 'x + ' num2str(P(2)) '  R^2 = ' num2str(rsq) '\n']);
%title(['Velocity and Vessel Diameter Dependancy | ' 'Fit: y = ' num2str(P(1)) 'x + ' num2str(P(2)) '  R^2 = ' num2str(rsq)])
title(['Velocity and Vessel Diameter Dependancy'])
figure(24)
clf
loglog(d,v,'Color',[0 0 1],'LineWidth',5);
hold on
loglog(exp(X),exp(Y),'k.','MarkerSize',1);
grid on
xlabel('Diameter (um)');
ylabel('Velocity (mm/s)');
title(['Velocity and Vessel Diameter Dependancy | ' 'Fit: y = ' num2str(P(1)) 'x + ' num2str(P(2)) '  R^2 = ' num2str(rsq)])
set(gca,'FontSize',12)
axis([5 45 0.01 10]);
yticklabels([0.01 0.1 1 10]);

figure(25)
clf
plot(log(d),log(v),'Color',[0 0 0],'LineWidth',5);
hold on
plot(X,Y,'Color',[0.1 0.1 1],'MarkerSize',1);
grid on
xlabel('Diameter (um)');
ylabel('Velocity (mm/s)');
title(['Velocity and Vessel Diameter Dependancy | ' 'Fit: y = ' num2str(P(1)) 'x + ' num2str(P(2)) '  R^2 = ' num2str(rsq)])
set(gca,'FontSize',12)
% axis([5 45 0.01 10]);
% yticklabels([0.01 0.1 1 10]);
%loglog(d,v_mean,'Color',[0.8 0 0],'LineWidth',5); % mean of data
%% Plot N dependancy with the vessel diameter
figure(96)
clf
hold on
longueur = 100;
for jj = 1:n_bubbles
    if(~isempty(bubbles{jj}.XYZ_laminar))
        if(length(bubbles{jj}.radii)>= longueur)
            title(jj);
            h = histogram(bubbles{jj}.radii(1:end));
            pause
        end
    end 
end
xlabel('d (um)');
ylabel('N');
%% Plot mean velocity and diameter dependancy
figure(22);
clf
set(gcf,'color','w');
colormap jet
plot(d,v,'LineWidth',5);
hold on
average = -1;
clear v_plot d_plot
for jj = 1:n_bubbles
    if(~isempty(bubbles{jj}.XYZ_laminar))
        bubbles{jj}.distances = sqrt(sum(diff(bubbles{jj}.XYZ_laminar,[],1).^2,2));
        v_plot(:,jj) = bubbles{jj}.distances./dt./1000;
        d_plot(:,jj) = 2*r_orig_smooth(bubbles{jj}.closest_nodes(2:end));
        hold on
    end 
end
axis([0 45 0 3.5]);
% xlabel('d (um)');
% ylabel('v (mm/s)');
% legend('Hingot','Gen data mean');
%% Plot radial distribution
figure(14)
clf
%scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled') % Shortest path nodes);
m = 2; % 
N = 0;
subplot(1,2,1)
for jj = 1:n_bubbles
    hold on
    if(size(bubbles{jj}.XYZ_laminar,1)>=m)
        plot1 = plot3(bubbles{jj}.XYZ_laminar(m,1),bubbles{jj}.XYZ_laminar(m,2),bubbles{jj}.XYZ_laminar(m,3),'.','Color', [(bubbles{jj}.poiseuille), 0, 1-bubbles{jj}.poiseuille]);
        plot1.Color(4) = 0.4;
        N = N+1;
    end
end
title('Radial distribution');
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
view(-95,-4)
subplot(1,2,2)
for jj = 1:n_bubbles
hold on
    if(size(bubbles{jj}.XYZ_laminar,1)>=m)
    plot1 = plot3(bubbles{jj}.XYZ_laminar(2,1),bubbles{jj}.XYZ_laminar(2,2),bubbles{jj}.XYZ_laminar(2,3),'.','Color', [(bubbles{jj}.poiseuille), 0, 1-bubbles{jj}.poiseuille]);
    plot1.Color(4) = 0.4;
    end
end
title('Radial distribution');
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
view(-100,90)
%% Plot radial distribution scatter3
figure(14)
clf
%scatter3(pos(:,1),pos(:,2),pos(:,3),1,[0 0 0],'filled') % Shortest path nodes);
m = 2;
n_bubbles = length(bubbles);
subplot(1,2,1)
for jj = 1:n_bubbles
    hold on
    bubbles{jj}.distances = sqrt(sum(diff(bubbles{jj}.XYZ_laminar,[],1).^2,2));
    if(size(bubbles{jj}.XYZ_laminar,1)>=m)
        n = length(bubbles{jj}.distances);
        scatter3(bubbles{jj}.XYZ_laminar(2,1), bubbles{jj}.XYZ_laminar(2,2), bubbles{jj}.XYZ_laminar(2,3),2,[bubbles{jj}.distances(2)],'filled')
        hold on
    end
end
title('Radial distribution');
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
view(-101,-2)
subplot(1,2,2)
for jj = 1:n_bubbles
hold on
    if(size(bubbles{jj}.XYZ_laminar,1)>=m)
        n = length(bubbles{jj}.distances);
        scatter3(bubbles{jj}.XYZ_laminar(2,1), bubbles{jj}.XYZ_laminar(2,2), bubbles{jj}.XYZ_laminar(2,3),1,[bubbles{jj}.distances(2)],'filled')
        hold on
    end
end
title('Radial distribution');
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
view(-100,90)

%% Poiseuille Validation
% min = 1;
% max = 0;
% idx_min = -1;
% idx_max = -1;
% for jj = 1:n_laminar
%     if(bubbles{jj}.poiseuille<min)
%         min = bubbles{jj}.poiseuille;
%         idx_min = jj;
%     elseif (bubbles{jj}.poiseuille>max)
%         max = bubbles{jj}.poiseuille;
%         idx_max = jj;
%     end
% end
% figure(13)
% clf
% for ii = 1:2
% hold on
%     % Fast bubble
%     jj = idx_min;
%     plot3(bubbles{jj}.XYZ_centerLine(1,:),bubbles{jj}.XYZ_centerLine(2,:),bubbles{jj}.XYZ_centerLine(3,:),'k', 'LineWidth',2) % Shortest path nodes);    
%     plot3(bubbles{jj}.XYZ_laminar(:,1),bubbles{jj}.XYZ_laminar(:,2),bubbles{jj}.XYZ_laminar(:,3),'o','Color', [(1-bubbles{jj}.poiseuille), 0, bubbles{jj}.poiseuille]);
%     jj = idx_max;
%     plot3(bubbles{jj}.XYZ_centerLine(1,:),bubbles{jj}.XYZ_centerLine(2,:),bubbles{jj}.XYZ_centerLine(3,:),'k', 'LineWidth',2) % Shortest path nodes);    
%     plot3(bubbles{jj}.XYZ_laminar(:,1),bubbles{jj}.XYZ_laminar(:,2),bubbles{jj}.XYZ_laminar(:,3),'o','Color', [(1-bubbles{jj}.poiseuille), 0, bubbles{jj}.poiseuille]);
% end
% legend('Ligne centrale','Bulle la plus rapide','Ligne centrale','Bulle la plus lente','location','best');
% xlabel('x','FontSize',20);
% ylabel('y','FontSize',20);
% zlabel('z','FontSize',20);
% title('Validation de l''algorithme de trajectorie laminaire');
%% Simulation of 1 trajectory
% v = 1000; %(um/s)
% dt = 0.01; % (s)
% inter_distance = v*dt; %(um)
% clear X Y Z points new_distances dd distances_point_previous distances_next_previous closest_nodes
% tic
% random_end_node = randi([1 length(end_nodes)],1,1);
% start = 2;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
% DG = digraph(s,t,r_inverse); % Directed graph generation
% [SP, ~] = shortestpathtree(DG,start,end_nodes(random_end_node)); % Shortest path
% edges = table2array(SP.Edges);
% min_bubbles = 10;
% if(length(edges)-4 > min_bubbles)
%     nodes = [edges(:,1);edges(end,2)]-1; % It is the previous node!
%     trajectory = pos(nodes,:); % Nodes positions attribution
%     distances = sqrt(sum(diff(trajectory,[],1).^2,2));
%     distances_cum = cumsum(distances);
%     d_trajectory = sum(sqrt(sum(diff(trajectory,[],1).^2,2)));
%     xyz = trajectory';
%     spline_f = cscvn(xyz);
%     coefficients = spline_f.coefs;
%     start = 0; % (um)
%     %new_distances = start:inter_distance:floor(d_trajectory);
%     %%% Vary the distance according to the closest node's radius
%     dd = 0;
%     new_distances(1,1) = start;
%     %previous_node_idxes(1,1) = 1;
%     %distances_point_previous(1,1) = start;
%     %distances_next_previous(1,1) = distances_cum(1)-start;
%     k = 2;
%     while(d_trajectory-new_distances(k-1,1) > inter_distance)
%         previous_nodes_idx = find(distances_cum <= dd); 
%         if(~isempty(previous_nodes_idx)) % if the starting distance is greater than the first node
%             previous_node_idxes(k-1,1) = previous_nodes_idx(end)+1;
%             distances_point_previous(k-1,1) = dd - distances_cum(previous_node_idxes(k-1,1)-1);
%             distances_next_previous(k-1,1) = distances_cum(previous_node_idxes(k-1,1)) - distances_cum(previous_node_idxes(k-1,1)-1);
%             if((distances_point_previous(k-1)/distances_next_previous(k-1))<=0.5) % if previous point is closest
%                 if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1))<d_trajectory) % if length not exceeding path length
%                     closest_nodes(k-1,1) = previous_node_idxes(k-1,1);
%                     new_distances(k,1) = dd + inter_distance*r_norm(closest_nodes(k-1));
%                     k = k+1;
%                 end
%             else
%                 if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1)+1)<d_trajectory) % if length not exceeding path length
%                     closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
%                     new_distances(k,1) = dd + inter_distance*r_norm(closest_nodes(k-1));
%                     k = k+1;
%                 end
%             end
%         else
%             previous_node_idxes(k-1,1) = 1;
%             next_node_idx = previous_node_idxes(k-1,1) + 1;
%             distances_point_previous(k-1,1) = dd;
%             distances_next_previous(k-1,1) = distances_cum(previous_node_idxes(k-1,1));
%             if((distances_point_previous(k-1)/distances_next_previous(k-1))<=0.5) % if previous point is closest
%                 if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1))<d_trajectory) % if length not exceeding path length
%                     closest_nodes(k-1,1) = previous_node_idxes(k-1,1);
%                     new_distances(k,1) = dd + inter_distance*r_norm(closest_nodes(k-1));
%                     k = k+1;
%                 end
%             else
%                 if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1)+1)<d_trajectory) % if length not exceeding path length
%                     closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
%                     new_distances(k,1) = dd + inter_distance*r_norm(closest_nodes(k-1));
%                     k = k+1;
%                 end
%             end
%         end
%         dd = new_distances(k-1,1);
%     end
%     L = length(new_distances)-1;
%     X = zeros(1,L);
%     Y = zeros(1,L);
%     Z = zeros(1,L);
%     points = zeros(3,L);
%     for ii = 1:L
%         d = new_distances(ii);
%         delta = distances_point_previous(ii)/sqrt(distances_next_previous(ii));%distance_point_previous/distance_next_previous; % distance_point_previous_normalized
%         d_euclidian = sqrt((trajectory(previous_node_idxes(ii)+1,1)-trajectory(previous_node_idxes(ii),1))^2+(trajectory(previous_node_idxes(ii)+1,3)-trajectory(previous_node_idxes(ii),2))^2+(trajectory(previous_node_idxes(ii)+1,3)-trajectory(previous_node_idxes(ii),3))^2);
%         % point calculation using spline
%         pp = (previous_node_idxes(ii)-1)*3 +1;
%         ax = coefficients(pp,1);
%         bx = coefficients(pp,2);
%         cx = coefficients(pp,3);
%         dx = coefficients(pp,4);
%         %delta = mean_segment_length+1;%(trajectory(previous_node_idx+1,1)-trajectory(previous_node_idx,1));
%         X(ii) = ax.*(delta.^3) + bx.*(delta.^2) + cx.*(delta) + dx;
%         ay = coefficients(pp+1,1);
%         by = coefficients(pp+1,2);
%         cy = coefficients(pp+1,3);
%         dy = coefficients(pp+1,4);
%         %delta = mean_segment_length+1;%(trajectory(previous_node_idx+1,2)-trajectory(previous_node_idx,2));
%         Y(ii) = ay.*(delta.^3) + by.*(delta.^2) + cy.*(delta) + dy;
%         az = coefficients(pp+2,1);
%         bz = coefficients(pp+2,2);
%         cz = coefficients(pp+2,3);
%         dz = coefficients(pp+2,4);
%         %delta = mean_segment_length+1;%(trajectory(previous_node_idx+1,3)-trajectory(previous_node_idx,3));
%         Z(ii) = az.*(delta.^3) + bz.*(delta.^2) + cz.*(delta) + dz;
%         points = [X;Y;Z];
%     end
%     toc
%     
% end     
% 
% figure(11)
% clf
% plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3),'.k') % Shortest path nodes);
% hold on
% [xyz_spline, values] = fnplt(spline_f);
% plot3(xyz_spline(1,:),xyz_spline(2,:),xyz_spline(3,:),'k');
% hold on
% plot3(points(1,:),points(2,:),points(3,:),'or');
% titre1 = '# originaux : ';
% titre2 = num2str(length(trajectory));
% titre3 = '| # générés : ';
% titre4 = num2str(length(points));
% titre_final = [titre1 titre2 ' ' titre3 titre4];
% title(titre_final);
% xlabel('x','FontSize',10);
% ylabel('y','FontSize',10);
% zlabel('z','FontSize',10);
% view(20,30)
% legend('Noeuds originaux','Lier les noeuds originaux','Points générés','location','best');
% set(gcf,'color','w');
%% Test trajectoires laminaires
% n_bubbles = 100;
% tic
% for jj = 1:n_bubbles
% clear xyz parallel perpendicular rayons
% xyz = points';
% parallel = [xyz(2:end,1)-xyz(1:end-1,1) xyz(2:end,2)-xyz(1:end-1,2) xyz(2:end,3)-xyz(1:end-1,3)];
% perpendicular = zeros(length(parallel),3);
% perpendicular2 = zeros(length(parallel),3);
% for i = 1:length(parallel)
%     perpendiculars = null(parallel(i,:)); 
%     perpendicular(i,:) = perpendiculars(:,1)'; %perpendicular vector 1
%     perpendicular2(i,:) = perpendiculars(:,2)'; %perpendicular vector 2
% end
% lin_combination = (-1+2*rand(1))*perpendicular+(-1+2*rand(1))*perpendicular2;
% % v1 = perpendicular(1,:);
% % v2 = perpendicular2(1,:);
% % d = 180;
% % c = 0;
% % theta = (d-c).*rand(1) + c;
% % x1 = cosd(theta)/norm(v1);
% % x2 = sind(theta)/norm(v2);
% % perpendicular(1,:) = norm(v1)*(x1*v1 + x2*v2); % Rotated vector
% 
% radii = r_orig_smooth(closest_nodes);
% compensation_radii = radii(2:end,1)-radii(1:end-1,1);
% compensation_rayon_cumul = cumsum(compensation_radii);
% %vR = perpendicular(:,1)*rand(1)*rayons(1);%% Ne pas oublier de changer ce vecteur en le rotatant
% laminar_xyz = zeros(length(xyz),3);
% laminar_xyz(1,:) = xyz(1,:)+ lin_combination(1,:).*radii(1);
% laminar_xyz(2:end,:) = xyz(2:end,:) + lin_combination.*radii(2:end);
% %laminar_xyz = xyz + ones(length(xyz),3).*perpendicular(1,:)*rayons(1);
% %laminar_xyz(2:end,:) = laminar_xyz(2:end,:)+ perpendicular.*compensation_rayon_cumulee;
%     laminar{jj}.trajectory = laminar_xyz;
% end
% elapsed = toc;
% figure(13)
% clf
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k','LineWidth',5) % Shortest path nodes);
% for jj = 1:n_bubbles
% hold on
%     plot1 = plot3(laminar{jj}.trajectory(:,1),laminar{jj}.trajectory(:,2),laminar{jj}.trajectory(:,3));
%     plot1.Color(4) = 0.4;
% end
% titre1 = 'Laminar flow simulation with ';
% titre2 = num2str(n_bubbles);
% titre3 = ' trajectories. Elapsed time : ';
% titre4 = num2str(elapsed);
% titre5 = ' s.';
% titre_final = [titre1 titre2 ' ' titre3 titre4 titre5];
% title(titre_final);
% legend('Center Line','Generated Laminar Flow','location','best');
% xlabel('x','FontSize',20);
% ylabel('y','FontSize',20);
% zlabel('z','FontSize',20);
% view(-75,20)
%% Vieille simulation
% clear bubbles trajectories DG
% %rng(0,'twister'); % initialize the random number generator
% v = 1000; %(um/s)
% dt = 0.01; % (s)
% inter_distance = v*dt; %(um)
% n_bubbles = 100;
% bubbles = cell(n_bubbles,1);
% finish = start + randi([1 size(s,1)-start-1],1,n_bubbles); % random integer from [start:#target_nodes-start]
% DG = digraph(s,t,r_inverse); % Directed graph generation
% param.radius_padding = 2.5; % To account for the size of the bubble
% v_max = 0;
% tot_toc = 0;
% min_bubbles = 100;
% jj = 1;
% while jj <= n_bubbles
%     clear X Y Z points xyz distances distances_cum new_distances dd distances_point_previous distances_next_previous
%     tic
%     random_end_node = randi([1 length(end_nodes)],1,1);
%     start = 2;%randi([1 size(s,1)],1,1) % random integer from [1:#source_nodes]
%     [SP, ~] = shortestpathtree(DG,start,end_nodes(random_end_node)); % Shortest path
%     edges = table2array(SP.Edges);
%     if(length(edges)-4 > min_bubbles)
%         bubbles{jj}.nodes = [edges(:,1);edges(end,2)]-1; % It is the previous node!
%         bubbles{jj}.trajectory = pos(bubbles{jj}.nodes,:); % Nodes positions attribution
%         %d_trajectory = sum(sqrt(sum(diff(bubbles{jj}.trajectory,[],1).^2,2))); % Path length
%         n_points = size(bubbles{jj}.trajectory,1);%floor(d_trajectory / d);
%         %n_points = round(d_trajectory/d) + 1;
%         %bubbles{jj}.spline_f = cscvn(bubbles{jj}.trajectory'); % Spline generation
%         %bubbles{jj}.trajectory_random = spline2Points(bubbles{jj}.spline_f.coefs,1); % Random points generation from spline coefficients
%         %bubbles{jj}.trajectory_random = interparc(n_points,bubbles{jj}.trajectory(:,1),bubbles{jj}.trajectory(:,2),bubbles{jj}.trajectory(:,3),'spline');
%         %%%%% 
%         distances = sqrt(sum(diff(bubbles{jj}.trajectory,[],1).^2,2));
%         distances_cum = cumsum(distances);
%         d_trajectory = sum(sqrt(sum(diff(bubbles{jj}.trajectory,[],1).^2,2)));
%         xyz = bubbles{jj}.trajectory';
%         spline_f = cscvn(xyz);
%         coefficients = spline_f.coefs;
%         start = 0; % (um)
%         %new_distances = start:inter_distance:floor(d_trajectory);
%         %%% Vary the distance according to the closest node's radius
%         dd = 0;
%         new_distances(1,1) = start;
%         %previous_node_idxes(1,1) = 1;
%         %distances_point_previous(1,1) = start;
%         %distances_next_previous(1,1) = distances_cum(1)-start;
%         k = 2;
%         while(d_trajectory-new_distances(k-1,1) > inter_distance)
%             previous_nodes_idx = find(distances_cum <= dd); 
%             if(~isempty(previous_nodes_idx)) % if the starting distance is greater than the first node
%                 previous_node_idxes(k-1,1) = previous_nodes_idx(end)+1;
%                 distances_point_previous(k-1,1) = dd - distances_cum(previous_node_idxes(k-1,1)-1);
%                 distances_next_previous(k-1,1) = distances_cum(previous_node_idxes(k-1,1)) - distances_cum(previous_node_idxes(k-1,1)-1);
%                 if((distances_point_previous(k-1)/distances_next_previous(k-1))<=0.5) % if previous point is closest
%                     if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1))<d_trajectory) % if length not exceeding path length
%                         bubbles{jj}.closest_nodes(k-1,1) = previous_node_idxes(k-1,1);
%                         new_distances(k,1) = dd + inter_distance*r_norm(bubbles{jj}.closest_nodes(k-1));
%                         k = k+1;
%                     end
%                 else
%                     if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1)+1)<d_trajectory) % if length not exceeding path length
%                         bubbles{jj}.closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
%                         new_distances(k,1) = dd + inter_distance*r_norm(bubbles{jj}.closest_nodes(k-1));
%                         k = k+1;
%                     end
%                 end
%             else
%                 previous_node_idxes(k-1,1) = 1;
%                 next_node_idx = previous_node_idxes(k-1,1) + 1;
%                 distances_point_previous(k-1,1) = dd;
%                 distances_next_previous(k-1,1) = distances_cum(previous_node_idxes(k-1,1));
%                 if((distances_point_previous(k-1)/distances_next_previous(k-1))<=0.5) % if previous point is closest
%                     if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1))<d_trajectory) % if length not exceeding path length
%                         bubbles{jj}.closest_nodes(k-1,1) = previous_node_idxes(k-1,1);
%                         new_distances(k,1) = dd + inter_distance*r_norm(bubbles{jj}.closest_nodes(k-1));
%                         k = k+1;
%                     end
%                 else
%                     if(dd + inter_distance*r_norm(previous_node_idxes(k-1,1)+1)<d_trajectory) % if length not exceeding path length
%                         bubbles{jj}.closest_nodes(k-1,1) = previous_node_idxes(k-1,1)+1;
%                         new_distances(k,1) = dd + inter_distance*r_norm(bubbles{jj}.closest_nodes(k-1));
%                         k = k+1;
%                     end
%                 end
%             end
%             dd = new_distances(k-1,1);
%         end
%         L = length(new_distances)-1;
%         X = zeros(1,L);
%         Y = zeros(1,L);
%         Z = zeros(1,L);
%         points = zeros(3,L);
%         for ii = 1:L
%             d = new_distances(ii);
%             delta = distances_point_previous(ii)/sqrt(distances_next_previous(ii));%distance_point_previous/distance_next_previous; % distance_point_previous_normalized
%             d_euclidienne = sqrt((bubbles{jj}.trajectory(previous_node_idxes(ii)+1,1)-bubbles{jj}.trajectory(previous_node_idxes(ii),1))^2+(bubbles{jj}.trajectory(previous_node_idxes(ii)+1,3)-bubbles{jj}.trajectory(previous_node_idxes(ii),2))^2+(bubbles{jj}.trajectory(previous_node_idxes(ii)+1,3)-bubbles{jj}.trajectory(previous_node_idxes(ii),3))^2);
%             % point calculation using spline
%             pp = (previous_node_idxes(ii)-1)*3 +1;
%             ax = coefficients(pp,1);
%             bx = coefficients(pp,2);
%             cx = coefficients(pp,3);
%             dx = coefficients(pp,4);
%             %delta = mean_segment_length+1;%(trajectory(previous_node_idx+1,1)-trajectory(previous_node_idx,1));
%             X(ii) = ax.*(delta.^3) + bx.*(delta.^2) + cx.*(delta) + dx;
%             ay = coefficients(pp+1,1);
%             by = coefficients(pp+1,2);
%             cy = coefficients(pp+1,3);
%             dy = coefficients(pp+1,4);
%             %delta = mean_segment_length+1;%(trajectory(previous_node_idx+1,2)-trajectory(previous_node_idx,2));
%             Y(ii) = ay.*(delta.^3) + by.*(delta.^2) + cy.*(delta) + dy;
%             az = coefficients(pp+2,1);
%             bz = coefficients(pp+2,2);
%             cz = coefficients(pp+2,3);
%             dz = coefficients(pp+2,4);
%             %delta = mean_segment_length+1;%(trajectory(previous_node_idx+1,3)-trajectory(previous_node_idx,3));
%             Z(ii) = az.*(delta.^3) + bz.*(delta.^2) + cz.*(delta) + dz;
%             bubbles{jj}.trajectory_random = [X;Y;Z]';
%         end
%         toc
%         %%%%%
%         bubbles{jj}.radius = r(bubbles{jj}.closest_nodes);   
%         bubbles{jj}.t = (0:dt:dt*(n_points-1))';
%         %****Points generation normal to spline ***
%         n = size(bubbles{jj}.trajectory_random,1);
%         bubbles{jj}.gen_points = zeros(n-1,3);
%         for i = 1:n-1
%             p2 = bubbles{jj}.trajectory_random(i+1,:);
%             p1 = bubbles{jj}.trajectory_random(i,:);
%             r1 = p2-p1; % Distances
%             r1 = r1/norm(r1); % Normalized distances
%             normal = null(r1(:).'); % Perpendicular vectors
%             v1 = normal(:,1)'; % Perpendicular vector 1
%             v2 = normal(:,2)'; % Perpendicular vector 2
%             c = 0; % Min angle
%             d = 180; % Max angle
%             a = -abs((bubbles{jj}.radius(i)-param.radius_padding)); % Min perpendicular amplitude
%             b = abs(bubbles{jj}.radius(i)-param.radius_padding); % Max perpendicular amplitude
%             theta = (d-c).*rand(1) + c;
%             x1 = cosd(theta)/norm(v1);
%             x2 = sind(theta)/norm(v2);
%             vR = norm(v1)*(x1*v1 + x2*v2); % Rotated vector
%             bubbles{jj}.r_bubble(i,1) = (b-a).*rand(1) + a; % Random radius given to the bubble
%             %bubbles{jj}.normalized_velocity(i) = 1 - bubbles{jj}.r_bubble(i).^2 ./ bubbles{jj}.radius(i).^2; % Normalized velocity following Poiseuille's law
%             bubbles{jj}.gen_points(i,:) = p1 + bubbles{jj}.r_bubble(i)*(vR); % Generated points normal to the trajectory
% 
%         end
%         bubbles{jj}.distances = sqrt(sum(diff(bubbles{jj}.gen_points,[],1).^2,2));
%         if(~isempty(bubbles{jj}.distances)) %if there are generated points
%             % New spline apposition to smooth data
%             xyz = (bubbles{jj}.gen_points)';
%             [ndim,npts]=size(xyz);
%             bubbles{jj}.smooth_xyz=zeros(size(xyz));
%             for k=1:ndim
%                bubbles{jj}.smooth_xyz(k,:)=ppval(csaps(1:npts,xyz(k,:),.5),1:npts);
%             end
%             % Velocity calculation
%             bubbles{jj}.velocities(1,1) = 0;
%             n_velocities = length(bubbles{jj}.gen_points);
%             bubbles{jj}.velocities(2:n_velocities) = (ones(n_velocities-1,1) - ((bubbles{jj}.r_bubble(1:n_velocities-1).^2)./((bubbles{jj}.radius(1:n_velocities-1)-param.radius_padding).^2))).* bubbles{jj}.distances ./ dt; % v = (R^2-r^2)* dx/dt
%         else
%             bubbles{jj}.velocities(1,1) = NaN;
%         end
%         v_local = max(bubbles{jj}.velocities);
%         if(v_local>v_max);v_max = v_local;end
%         bubbles{jj}.velocities_normalized = bubbles{jj}.velocities./v_max;
%         progress = jj/n_bubbles*100;
%         %fprintf('Le processus est à %6.5f%s\n',progress,'%');
%         jj = jj+1;
%     %jj
%     else
%         % do nothing    
%     end
%     %tot_toc = DisplayEstimatedTimeOfLoop(tot_toc+toc, jj, n_bubbles);
% end
% %save('bubbles.mat','bubbles','-v7.3');
% 
% % Affichage
% tic
% figure(9);
% clf
% hold on
% %plot3(pos(:,1),pos(:,2),pos(:,3),'.k')
% view(10,20)
% n_bubbles = size(bubbles,1);
% for jj = 1:n_bubbles
%     %plot3(bubbles{jj}.interpolation(:,1),bubbles{jj}.interpolation(:,2),bubbles{jj}.interpolation(:,3),'.r')
%     %plot3(bubbles{jj}.trajectory(:,1),bubbles{jj}.trajectory(:,2),bubbles{jj}.trajectory(:,3))
%     %plot3(bubbles{jj}.trajectory_random(:,1),bubbles{jj}.trajectory_random(:,2),bubbles{jj}.trajectory_random(:,3),'.r')
%     %plot3(bubbles{jj}.gen_points(:,1), bubbles{jj}.gen_points(:,2), bubbles{jj}.gen_points(:,3),'.','MarkerSize',1);
%     if(~isempty(bubbles{jj}))
%         n = length(bubbles{jj}.velocities_normalized);
%         if(~isempty(bubbles{jj}.distances)&&size(bubbles{jj}.smooth_xyz,2)==size(bubbles{jj}.velocities_normalized,2))
%             scatter3(bubbles{jj}.smooth_xyz(1,:), bubbles{jj}.smooth_xyz(2,:), bubbles{jj}.smooth_xyz(3,:),1,[zeros(n,1) bubbles{jj}.velocities_normalized' ones(n,1) - bubbles{jj}.velocities_normalized']);
%         end
%     end
% end
% %plot3(bubbles{jj}.gen_points(:,1), bubbles{jj}.gen_points(:,2), bubbles{jj}.gen_points(:,3),'.','MarkerSize',1);
% xlabel('x','FontSize',20);
% ylabel('y','FontSize',20);
% zlabel('z','FontSize',20);
% titre = ['Simulation of ',num2str(n_bubbles),' bubbles - interparc'];
% title(titre);
% legend('Bubbles (Blue = slow) (Green = fast)','location','best')
% toc
%% Création d'une vidéo
% figure
% volshow(seg);
% OptionZ.FrameRate=15;OptionZ.Duration=10;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'V_Surface',OptionZ)
% V(1:50, 1:60, 1:800) = 255;

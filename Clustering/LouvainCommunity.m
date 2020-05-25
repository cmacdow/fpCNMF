function [ovr_comm, ovr_Q] = LouvainCommunity(w, varargin),

%Expects a weighted graph (not directional for now).  Should be NxN where
%w(i, j) is the weight from i to j.

%Define options
opts.MaxIterations = 25; %maximum # of iterations
opts.NumRandomStarts = 5; %number of random starts
opts.RandOrder = 1; %use random order?
opts.Verbose = 1; %report progress

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end

%% Check inputs

if ~opts.RandOrder & (opts.NumRandomStarts > 1),
    warning('Not worth doing random restarts without random index orders.');
    opts.NumRandomStarts = 1;
end

%% Loop through network, building communities

ovr_Q = 0;
ovr_comm = [];

%Network variables that won't change
two_m = sum(w(:));

for cur_rnd = 1:opts.NumRandomStarts,
    if opts.Verbose,
        fprintf('\nStarting clustering attempt %d/%d...\n', cur_rnd, opts.NumRandomStarts);
    end
    
    %Initialize cur_w
    cur_w = w;
    
    %Initialize the network clusters
    N = size(cur_w, 1);
    cur_comm = [1:N]'; %should be N x 1
    %Matrix for tracking clusters over iterations
    comm = cur_comm;
    
    %Calculate current modularity
    cur_Q = calcModularity(cur_w, cur_comm);
    
    m = sum(sum(cur_w));
    MOD = 0;
    COMu = unique(cur_comm);
    for j=1:length(COMu)
        Cj = find(cur_comm==COMu(j));
        Ec = sum(sum(cur_w(Cj,Cj)));
        EtIn = sum(sum(cur_w(Cj,:)));
        EtOut = sum(sum(cur_w(:,Cj)));
        MOD = MOD + Ec/m-EtIn*EtOut/m^2;
    end
    
    %Update the network recursively
    cur_iter = 0;
    stop = 0;
    while ~stop,
        %Count number of iterations
        cur_iter = cur_iter + 1;
        
        %Initialize communities for current network
        N = size(cur_w, 1);
        cur_comm = [1:N]';
        next_Q = cur_Q;
        
        if opts.Verbose,
            fprintf('\tIteration %d: %d nodes; modularity = %5.4f...\n', cur_iter, N, cur_Q);
        end
        
        %Initialize terms that are constant for a given network
        %Calculate weighted in and out degree of each node
        k_in = sum(cur_w, 1);
        k_out = sum(cur_w, 2);
        %Sum total inputs and outputs into a node (community)
        sum_tot_in = sum(cur_w, 1);
        sum_tot_out = sum(cur_w, 2);
        
        %Loop through, updating cluster identities until we can't improve modularity
        clust_stop = 0; clust_iter = 0;
        while ~clust_stop,
            %Set stop to true, will revert to false once change is made
            clust_stop = 1;
            
            %Generate list of nodes to visit
            if opts.RandOrder, 
                ind_list = randperm(N);
            else
                ind_list = [1:N];
            end
            
            %Re-lable communities to go from 1 to # communities
            [~, ~, cur_comm] = unique(cur_comm); %
            uniq_comm = unique(cur_comm);
            %Count communities
            num_clusters = length(uniq_comm);
            %Keep community list
            clust_ind = cell(num_clusters, 1);
            sum_tot_in = zeros(num_clusters, 1);
            sum_tot_out = zeros(num_clusters, 1);
            for cur_clust = 1:num_clusters,
                clust_ind{cur_clust} = find(cur_comm == uniq_comm(cur_clust))';
                sum_tot_in(cur_clust) = sum(k_in(clust_ind{cur_clust}));
                sum_tot_out(cur_clust) = sum(k_out(clust_ind{cur_clust}));
            end
            
            if opts.Verbose,
                fprintf('\t\tClustering iteration %d: %d clusters, modularity = %6.5f\n', clust_iter, num_clusters, next_Q);
            end
            
            %Loop through list, attempting to cluster nodes
            for ind = 1:N,
                i = ind_list(ind);
                
                %Remove this node from its community, set default
                sum_tot_in(cur_comm(i)) = sum_tot_in(cur_comm(i)) - k_in(i);
                sum_tot_out(cur_comm(i)) = sum_tot_out(cur_comm(i)) - k_out(i);
                clust_ind{cur_comm(i)} = setdiff(clust_ind{cur_comm(i)}, i);
                               
                %Find all nodes that receive a connection from or send a connection to this node
                j_list = (cur_w(i, :) > 0) | (cur_w(:, i) > 0)';
                %Find the communities those belong to
                test_comm = unique(cur_comm(j_list));
                
                %Calculate the improvement in modularity for each possible community
                delta_Q = zeros(length(test_comm), 1);
                for cur_clust_ind = 1:length(test_comm),
                    %cur_clust = find(uniq_comm == test_comm(cur_clust_ind), 1, 'first');
                    cur_clust = test_comm(cur_clust_ind);
                    
                    %Calculate weight between this node and nodes in putative new cluster
                    k_i_in = sum(cur_w(i, clust_ind{cur_clust})) + sum(cur_w(clust_ind{cur_clust}, i));
                    
                    % Calculate current delta Q; modified for directed graphs
                    delta_Q(cur_clust_ind) = k_i_in - (k_in(i)*sum_tot_out(cur_clust) + k_out(i)*sum_tot_in(cur_clust))/two_m;
                end %loop of test communities
                delta_Q = delta_Q/two_m;
                
                
                %Test whether modularity increased with any of these moves
                [max_delta_Q, max_Q_ind] = max(delta_Q);
                if (max_delta_Q > 0) & (test_comm(max_Q_ind) ~= cur_comm(i)),
                    %Move node i to new community, need to update parameters for
                    %this community
                    prev_comm = cur_comm(i);
                    cur_comm(i) = test_comm(max_Q_ind);
                    if opts.Verbose & (mod(i-1, floor(N/3)) == 0),
                        fprintf('\t\t\tMoving %d from community %d to %d for %2.1e Q.\n', i, prev_comm, cur_comm(i), max_delta_Q);
                    end
                    
                    %Update modularity
                    k_i_in = sum(cur_w(i, clust_ind{prev_comm})) + sum(cur_w(clust_ind{prev_comm}, i));
                    prev_delta_Q = k_i_in - (k_in(i)*sum_tot_out(prev_comm) + k_out(i)*sum_tot_in(prev_comm))/two_m;  
                    prev_delta_Q = prev_delta_Q/two_m;
                    next_Q = next_Q + max_delta_Q - prev_delta_Q;
                    
                    %We made a change, so don't stop with this iteration
                    clust_stop = 0;
                end
                
                %Update our cluster statistics
                sum_tot_in(cur_comm(i)) = sum_tot_in(cur_comm(i)) + k_in(i);
                sum_tot_out(cur_comm(i)) = sum_tot_out(cur_comm(i)) + k_out(i);
                clust_ind{cur_comm(i)} = cat(2, clust_ind{cur_comm(i)}, i);
                
            end %node loop
            clust_iter = clust_iter + 1;
        end %cluster loop
        
        
        %% Build new network
        if opts.Verbose,
            fprintf('\tBuilding a new network.\n');
        end
            
        %Combine connections within a community into a single node
        N = num_clusters;
        new_w = sparse(N, N);
        for i = 1:N,
            for j = 1:N,
                new_w(i, j) = sum(sum(cur_w(clust_ind{i}, clust_ind{j})));
            end
        end
        
        % Initialize new communities
        new_comm = [1:N]';
        
        %Calculate new modularity
        new_Q = calcModularity(new_w, new_comm);       
        
        %% Check if modularity has improved
        if new_Q > cur_Q,
            %Convert the community to original nodes
            next_comm = cur_comm(comm(:, 1));
            
            % Copy to counter variables
            comm = cat(2, next_comm, comm);
            cur_Q = new_Q;
            cur_w = new_w;
            cur_comm = new_comm;
            
        elseif (cur_Q <= new_Q),
            %Time to stop, modularity isn't getting any better
            stop = 1;
        end
        
        %Check if we've run maximum number of iterations
        if (cur_iter >= opts.MaxIterations),
            stop = 1;
        end
        
    end %iteration loop
    
    %Check to see if this is better than overall best
    if (cur_Q > ovr_Q),
        ovr_comm = comm(:, 1:(end-1)); %take out the identity clustering
        ovr_Q = cur_Q;
    end
    
end %random loops
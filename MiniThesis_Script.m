%%% minithesis script

%% This function return the matrix after thresholding

function [matrix]= Thresholding(matrix)

	% Because some partial correlations from the synthetic data exceed the range [-1, 1], I turn those values into 1
	matrix(matrix> 1)= 1;
	% Because some partial correlations values are poorly relevant, I put a threshold to 0.05
	matrix( matrix< 0)= 0;

end

%% This function returns the degree of each node of the matrix

function [degrees]= NodeDegree(matrix)
	degrees= sum(matrix, 2);
	
end

%% This function returns the betweeness of each node in the network

function [betweenness]= BetweennessCentrality(matrix)
	% computing the characteristic pathlenght
	matrix(matrix> 0)= 1./matrix(matrix> 0); % inverting the entries, so the hiegher correlation nodes will have the shorter pathlenght
	d= distance_wei(matrix); 
    
	% computing the betwenness
	betweenness= betweenness_wei(d);
end

%% This function returns the Closeness of each node of the network

function [closeness]= ClosenessCentrality(matrix)
	% computing the distance of each node
	matrix(matrix> 0)= 1./matrix(matrix> 0); % inverting the entries, so the hiegher correlation nodes will have the shorter pathlenght
	d= distance_wei(matrix);
	% computing the sum of distances 
	d_sum= sum(d, 2);
	
	% computing Closeness
	closeness= 1./d_sum;
end

%% This function returns the z-scores of each node, for intramodule connectivity

function [zscores]= ZScores(matrix)
	% computing the modularity of the matrix
	[M, Q]= modularity_und(matrix);
	% computing the z-scores
	zscores= module_degree_zscore(matrix, M);
end

%% This function compute the participation coefficient of each node

function [participation]= ParticipationCoefficient(matrix)
	% computing the modularity of the matrix
	[M, Q]= modularity_und(matrix);
	% computing the participation ocefficient
	participation= participation_coef(matrix,M);
end

%% This function evalate the matrix and output the hub-scores of all the nodes of the matrix in the same order in which they are presented on the matrix.

function [hubs, scores]= Hubs(matrix)

	scores= zeros(1, size(matrix,1));
	top_25 = round(0.25*size(matrix,1));

	% scores for the top 25% higher degree
	degrees= NodeDegree(matrix);
	[~, I]= sort(degrees, 'descend');
   
	scores(I(1:top_25))= scores(I(1:top_25)) + 1;
    

	% scores for the top 25% higher betweenness centrality
	betweenness= BetweennessCentrality(matrix);
	[~, I]= sort(betweenness, 'descend');
    
	scores(I(1:top_25))= scores(I(1:top_25)) + 1;

	% scores for the top 25% higher closenmess centrality
	closeness= ClosenessCentrality(matrix);
	[~, I]= sort(closeness, 'descend');
	scores(I(1:top_25))= scores(I(1:top_25)) + 1;

	% scores for the top 25% higher z-score
	zscores= ZScores(matrix);
	[~, I]= sort(zscores, 'descend');
	scores(I(1:top_25))= scores(I(1:top_25)) + 1;

	% scores for the top 25% higher participation coefficient binary
	participation= ParticipationCoefficient(matrix);
	[~, I]= sort(participation, 'descend');
	scores(I(1:top_25))= scores(I(1:top_25)) + 1;
	
	% finding hubs
	[~, I]= sort(scores, 'descend');
    hubs(:,1)= I(1:top_25);
	

end

%% This function returns the benchmark of differneces

function [differences]= PermutationTest(matrixA, matrixB, NumberOfPermutations)

    differences= zeros(NumberOfPermutations, 1);
    A= squareform(matrixA);
	IA= find(A> 0);
	n= numel(A(IA));

    names= ["FA-Fs"; "FA-Con"; "FA-Mot"; "FA-Act"; "Cogni"; "Depr"; "PhF";  "SocF"; "RolePh"; "RoleE"; "EmotWB"; "Pain"; "HealthP"; "ChangeH"; "FutU"; "Visual";  "Motor"; "CommD";  "HA"; "Seiz"; "Drow"];

	
	B= squareform(matrixB);
    IB= find(B> 0);
    m= numel(B(IB));
	
	
     values = [A(IA), B(IB)];

    for i= 1:NumberOfPermutations

        newA = zeros(size(A));
        newB = zeros(size(B));

        I= randperm(n+m);
		
		newA(IA)= values(I(1:n));
		newA= squareform(newA);
		[hubsA, ~]= Hubs(newA);
        hubsA= names(hubsA);
		
		newB(IB)= values(I(n+1:end));
		newB= squareform(newB);
		[hubsB, ~]= Hubs(newB);	
    	hubsB= names(hubsB);

		
        differences(i,1)= 5- sum(ismember(hubsA, hubsB));
     end
 end

%% This function returns the p-value of significance of each hubs in a set of hubs

function [hubs]= Significance(matrix, NumberOfRandomizations)
    
    %computing the hubs and hub-scores of the matrix
    [hubs, OriginalScores]= Hubs(matrix);

    %computing the null-benchmark of randomized networks
    count = zeros(size(matrix, 1), 1);
    ScoresFromRandomizaiton= zeros(size(matrix,1), NumberOfRandomizations);
    A= squareform(matrix);
	IA= find(A> 0);
	n= numel(A(IA));
    values= [A(IA)];

    for i= 1:NumberOfRandomizations

        newA = zeros(size(A));
       
        I= randperm(n);
		
		newA(IA)= values(I(1:end));
		newA= squareform(newA);
		[~, ScoresFromRandomizaiton(:,i)]= Hubs(newA);

        count = count + (ScoresFromRandomizaiton(:, i) >= OriginalScores');
        

    end

    % computing the p-value for each hub
    pvs= count./ NumberOfRandomizations;

    hubs(:,2)= pvs(hubs(:,1));


	names= ["FA-Fs"; "FA-Con"; "FA-Mot"; "FA-Act"; "Cogni"; "Depr"; "PhF";  "SocF"; "RolePh"; "RoleE"; "EmotWB"; "Pain"; "HealthP"; "ChangeH"; "FutU"; "Visual";  "Motor"; "CommD";  "HA"; "Seiz"; "Drow"];

    for i= 1:size(hubs,1)

        figure;
        histogram(nonzeros(ScoresFromRandomizaiton(hubs(i,1),:)));
        
        % Add title and labels
        title(['Randomized hub score distribution of ', names(hubs(i,1))]);
        xlabel('Hub Scores');
        ylabel('Frequency');
        
        % Add vertical line (e.g., at mean)
        xline(OriginalScores(hubs(i,1)), 'r--', 'LineWidth', 3);
        
        % Set x-axis limits
        xlim([0, 6]); 
        
    end




end



%% Start of the actual code
names= ["FA-Fs"; "FA-Con"; "FA-Mot"; "FA-Act"; "Cogni"; "Depr"; "PhF";  "SocF"; "RolePh"; "RoleE"; "EmotWB"; "Pain"; "HealthP"; "ChangeH"; "FutU"; "Visual";  "Motor"; "CommD";  "HA"; "Seiz"; "Drow"];
%importing the matrices
matrixA= table2array(readtable("n_lowgrade"));
matrixB= table2array(readtable("n_highgrade"));

% Thresholding the matrices
matrixA= Thresholding(matrixA);
matrixB= Thresholding(matrixB);

figure; imagesc(matrixA); colorbar;
set(gca, 'XTick', 1:length(names), 'XTickLabel', names);
set(gca, 'YTick', 1:length(names), 'YTickLabel', names);
xtickangle(45);
title('Heatmap of the non-fatigued patients'' matrix');

figure; imagesc(matrixB); colorbar;
set(gca, 'XTick', 1:length(names), 'XTickLabel', names);
set(gca, 'YTick', 1:length(names), 'YTickLabel', names);
xtickangle(45);
title('Heatmap of the fatigued patients'' matrix');

% Computing the significance of these hubs
hubsA= Significance(matrixA, 1000);

for i= 1:size(hubsA,1)
disp([names(hubsA(i,1)), 'is a hub of the non-fatigued patients'' matrix, with a p-value of ', num2str(hubsA(i,2))]);
end

hubsB= Significance(matrixB, 1000);

for i= 1:size(hubsB,1)
disp([names(hubsB(i,1)), 'is a hub of the fatigued patients'' matrix, with a p-value of ', num2str(hubsB(i,2))]);
end

% computing the hubs difference
originalDifference= 5- sum(ismember(names(hubsA(:,1)), names(hubsB(:,1))));

% computing permutation test
differences= PermutationTest(matrixA, matrixB, 1000);
count= sum(differences>= originalDifference);
pv= count/ 1000;

figure;
histogram(differences);

% Add title and labels
title('Distribution of randomized comparisons');
xlabel('Number of different hubs per comparison');
ylabel('Frequency');

% Add vertical line (e.g., at mean)
xline(originalDifference, 'r--', 'LineWidth', 3);

% Set x-axis limits
xlim([-1, 6]); 

disp(['The original difference in hubs has a p-value of ', num2str(hubsB(i,2))]);
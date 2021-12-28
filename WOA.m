function [leaderScore, leaderPos, convergenceCurve] = WOA(noSearchAgents, noUsers, noSubcs, maxIter, lb, ub, fobj, X) 
	leaderPos = zeros(noUsers, noSubcs); 
	leaderScore = inf; 

	leader_score_pre = leaderScore; 

	convergenceCurve = zeros(1, maxIter); 

	% ======================== Initialization =================
	positions = zeros(noUsers, noSubcs, noSearchAgents);

	% If each variable has a different lb and ub
	for i = 1:noSearchAgents
		positions(:, :, i) = rand(noUsers, noSubcs).*(ub - lb) + lb; 
    end
    
	% Loop counter 
	t = 0;
	todoTol = 1; 
	delta = 1e-5; 
	flag = 0; 

	% Main loop
	while t < maxIter && flag < 7 

		for i = 1:noSearchAgents
			% Return back the search agents that go beyond the boundaries of the search space
			tmp = positions(:, :, i); 
			flag4lb = tmp < lb; 
			flag4ub = tmp > ub; 
			positions(:, :, i) = tmp.*(~(flag4lb + flag4ub)) + lb.*flag4lb + ub.*flag4ub; 

			% Calculate objective function for each search agent
			fitness = fobj(positions(:, :, i), X); 

			% Update the leader 
			if fitness < leaderScore
				leaderScore = fitness; 
				leaderPos = positions(:, :, i);
			end 
		end

		% a decreases linearly from 2 to 0
		a = 2 - t*(2/maxIter); 		
		% a2 linearly decreases from -1 to -2 to calculate t 
		a2 = -1 + t*(-1/maxIter); 

		% Update the position of each search agents
		for k = 1:noSearchAgents
			r1 = rand(); 
			r2 = rand(); % do we actually need r2  
			
			A = 2*a*r1 - a; 
			C = 2*r2;

			% parameters for spiral updating position
			b = 1; 
			l = (a2 - 1)*rand + 1; 

			p = rand(); 

			for i = 1:noUsers
				for j = 1:noSubcs
					% follow the shrinking encircling mechanism or prey search
					if p < 0.5
						% search for prey (exploration phase)
						if abs(A) >= 1
							randLeaderIndex = floor(noSearchAgents*rand + 1); 
							X_rand = positions(:, :, randLeaderIndex); 
							D_X_rand = abs(C*X_rand(i, j) - positions(i, j, k)); 
							positions(i, j, k) = X_rand(i, j) - A*D_X_rand; 
						elseif abs(A) < 1
							D_Leader = abs(C*leaderPos(i, j) - positions(i, j, k)); 
							positions(i, j, k) = leaderPos(i, j) - A*D_Leader; 
						end
					elseif p >= 0.5
						distance2Leader = abs(leaderPos(i, j) - positions(i, j, k)); 
						positions(i, j, k) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + leaderPos(i, j); 
					end 
				end 
			end 
		end

		% increase the iteration index by 1 
		t = t + 1; 
		convergenceCurve(t) = leaderScore; 

		if todoTol == 1 && abs(leaderScore - leader_score_pre) < delta
			flag = flag + 1; 
			convergenceCurve = convergenceCurve(1, 1:t);
		else 
			flag = 0; 
        end 
        [t leaderScore flag]; 
		leader_score_pre = leaderScore;
	end


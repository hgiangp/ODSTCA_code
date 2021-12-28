%--------------------------
% BWOA combining with two conditions to solve ODSTCA pb (Algorithm 2)
%---------------------------

function [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, conver_curve, conver_curve_woa, no_WOA_run, time] = BWOA2(functionName, doTol, noSearchAgents, noUsers, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa, gamma, eta, hArray, n0, pmin, pmax, Adet)
	tic 
	% initialization ======
	positions = Initialization(functionName, noUsers, noSubcs, noSearchAgents, Adet);  

	leader_pos_bwoa = zeros(noUsers, noSubcs); 
	leader_score_bwoa = -inf; 
	leader_score_pre = leader_score_bwoa; 
	% leader_score_woa 
	leader_pos_woa = zeros(noUsers, noSubcs); 


	% loop counter
	todoTol = doTol; 
	delta = 1e-5; 
	flag = 0; 
	
	conver_curve = zeros(1, maxIter);
    conver_curve_woa = zeros(1, maxIter); 
	C_Step = zeros(noUsers, noSubcs); 
	iter = 0; 
	phi = @(y,a,x,eta) y*log2(1 + a*x) - (a/log(2))*(eta + y*x)/(1 + a*x);
    fmin = @(y,a,x,eta) (eta+y*x)/(log2(1 + a*x)); 
	

	varepsilon = 1e-5 ;
 
	no_WOA_run = 0; 

	while iter < maxIter
 %       display(['iter = ' num2str(iter)]);  
		for k = 1:noSearchAgents

            fitbwoa = fobj_bwoa(positions(:, :, k)); %NS
 %           display([k fitbwoa score_bwoa leader_score_bwoa]);
 
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% define constraint to cut down on the number of WOA runs
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			if fitbwoa > 1e2 || fitbwoa <= leader_score_bwoa
				continue
			end 

			WOA_tmp = 0; 
            p_tmp = zeros(noUsers, noSubcs); 
			for i = 1:noUsers
				for j = 1:noSubcs
					if positions(i, j, k) == 0
						continue; 
					end 
					if phi(gamma(i), hArray(i, j)/n0, pmax, eta(i)) <= 0
						p_tmp(i, j) = pmax; 
					else
						p_s = pmin; p_t = pmax; 
						while (abs(p_t - p_s) > varepsilon)
							p_l = (p_t + p_s)/2; 
							if phi(gamma(i), hArray(i, j)/n0, p_l, eta(i)) <= 0
								p_s = p_l;
							else 
								p_t = p_l; 
							end 
						end 
						p_tmp(i, j) = (p_s + p_t)/2; 
					end
					WOA_tmp = WOA_tmp + fmin(gamma(i), hArray(i, j)/n0, p_tmp(i, j), eta(i));  
				end
			end 


            % update 30/3/2020 
            if (fitbwoa - WOA_tmp) <= leader_score_bwoa
                continue
            end 

            no_WOA_run = no_WOA_run + 1; 

            [WOA_rs, pos_woa, ~] = WOA(noSearchAgents, ...
				noUsers, noSubcs, 150, lb_woa, ub_woa, fobj_woa, positions(:, :, k));

            fitness = fitbwoa - WOA_rs; 

			% update the leader
			if fitness > leader_score_bwoa
				leader_score_bwoa = fitness; 
				leader_pos_bwoa = positions(:, :, k); 
				leader_pos_woa = pos_woa; 
				leader_score_woa = WOA_rs; 
			end 
		end

		a = 2 - iter*(2/maxIter); % a decreases linearly from 2 to 0
		a2 = -1 + iter*(-1/maxIter); % a2 decreases linearly from -1 to -2 
		% update the position of each search agents 

		for k = 1:noSearchAgents
			r1 = rand(); 
			r2 = rand(); 
			A = 2*a*r1 - a; 
			C = 2*r2; 
			% parameters for spiral updating position
			b = 1; 
			l = (a2 - 1)*rand + 1; 
			p = rand();
			for i = 1:noUsers
				for j = 1:noSubcs
					if p < 0.5
						% search for prey (exploration phase)
						if abs(A) >= 1 
							rand_idx = floor(noSearchAgents*rand + 1); 
							X_rand = positions(:, :, rand_idx); 
							D_X_rand = abs(C*X_rand(i, j) - positions(i, j, k)); 
							C_Step(i, j) = X_rand(i, j) - A*D_X_rand;
						elseif abs(A) < 1
							% shrinking encircling mechanism (exploitation phase)
							D_leader = abs(C*leader_pos_bwoa(i, j) - positions(i, j, k)); 
							C_Step(i, j) = leader_pos_bwoa(i, j) - A* D_leader; 
						end
					elseif p >= 0.5
						distance2leader = abs(leader_pos_bwoa(i, j) - positions(i, j, k));
						C_Step(i, j) = distance2leader*exp(b.*l).*cos(l.*2*pi) + leader_pos_bwoa(i, j); 
					end 

					sigmoid = 1/(1 + exp(-10*(C_Step(i, j)-0.5))); 

					p_rand = rand(); 
					if p_rand < sigmoid
						positions(i, j, k) = ~positions(i, j, k); 
					end 		
				end 
			end 
		end 
		iter = iter + 1; 
		conver_curve(iter) = leader_score_bwoa; 
		conver_curve_woa(iter) = leader_score_woa;

		[iter leader_score_woa leader_score_bwoa]

		if todoTol == 1 && iter > 40 && abs(leader_score_bwoa - leader_score_pre) <= delta 
			flag = flag + 1; 
		else 
			flag = 0; 
		end 
		leader_score_pre = leader_score_bwoa; 
		if flag == 15
			conver_curve(iter+1:maxIter) = conver_curve(iter); 
			conver_curve_woa(iter+1:maxIter) = conver_curve_woa(iter); 
			break; 
        end
	end 
	toc; 
	time = toc; 
end 
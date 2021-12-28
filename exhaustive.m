%---------------------
% Algorithm to solve ODSTCA pb by
% using EXHAUSTIVE SEARCH to solve SA pb.  
%---------------------
function [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, time] = exhaustive(noSearchAgents, noUsers, noSubcs, maxIter, lb, ub, fobj, lambda, f_l, f0)
	tic 
	sa = zeros(noUsers, noSubcs); 
	i = 1; 
    cnt = 0; 
	  
	leader_pos_bwoa = zeros(noUsers, noSubcs); 
	leader_score_bwoa = -inf; 
	leader_pos_woa = zeros(noUsers, noSubcs); 

	TRY(i); 

	function [] = solution()
		% TRANSMIT POWER ALLOCATION 
		[woa_rs, woa_pos, ~] = WOA(noSearchAgents, noUsers, noSubcs, maxIter, lb, ub, fobj, sa); 
		% COMPUTING RESOURCE ALLOCATION 
		dec = sum(sa, 2); 
		tmp = sqrt(lambda(:, 1).*f_l); 
		fobj3 = sum(dec.*tmp)^2/f0;
		fobj1 = sum(sum(lambda, 2).*dec); 
		bwoa = fobj1 - woa_rs - fobj3; 
		cnt = cnt + 1;  
		if bwoa > leader_score_bwoa
			leader_score_bwoa = bwoa; 
			leader_pos_bwoa = sa; 
			leader_pos_woa = woa_pos; 
		end 

	end 

	function [] = TRY(i)
		for j = 0:noSubcs 
			if j ~= 0
				sa(i, j) = 1; 
			end
			if i == noUsers
				solution();
			else 
				TRY(i + 1);
			end
			if j ~= 0
				sa(i, j) = 0; 
			end 

		end 
    end
    toc;  
    time = toc; 
end 

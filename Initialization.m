% --------------------
% Intialization the position of whales (subchannel assignment) 
% --------------------

function [positions] = Initialization(functionname, noUsers, noSubcs, noSearchAgents, Adet)
	positions = zeros(noUsers, noSubcs, noSearchAgents); 
	switch functionname
		case 'ODSTCA'
			for k = 1:noSearchAgents
				for i = 1:noUsers
					if rand > 0.5
						rand_idx = floor(noSubcs*rand + 1); 
						positions(i, rand_idx, k) = 1; 
					end 
				end 
			end
		case 'ARJOA'
			for k = 1:noSearchAgents
				for i = 1:noUsers
					rand_idx = floor(noSubcs*rand + 1); 
					positions(i, rand_idx, k) = 1; 
				end 
			end
		case 'IOJOA'
			for k = 1:noSearchAgents
				for i = 1:noUsers
					if Adet(i) == 1
						rand_idx = floor(noSubcs*rand + 1); 
						positions(i, rand_idx, k) = 1;
					end 
				end 	 
			end 

		case 'OFDMA'
			for k = 1:noSearchAgents
				for j = 1:noSubcs
					if rand > 0.5
						rand_jdx = floor(noUsers*rand + 1); 
						positions(rand_jdx, j, k) = 1; 
					end 
				end 
			end 
	end
end 
	

%-----------------------------------
% Obtain the objective function and search domain for WOA and BWOA. 
%-----------------------------------

function [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails(F, noUsers, noSubcs, hArray, p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adet)
% function name: F: depends on each simulation schemes {ARJOA, ODSTCA, IOJOA, OFDMA} 
% network size: noUser x noSubs
% channel gain: hArray 
% users' power budget: [p_min, p_max]
% tolerant power for NOMA: P_tol
% thermal noise: n0
% eta, gamma, f_l: parameters of local computing 
% penalty factor for constraint dealing: mu
% perference in time and energy: lambda = [lambda_t lambda_e]: 
% matrix defines offloading decision of IOJOA: Adet

	lb_woa = zeros(noUsers, noSubcs); 
	ub_woa = zeros(noUsers, noSubcs);

	PMin = zeros(noUsers, 1); 
	PMax = zeros(noUsers, 1); 

	hArray = hArray; 

	PMin(:) = p_min; 
	PMax(:) = p_max; 

	switch F
	case {'ODSTCA'}
		fobj_woa = @FWOA; 
		lb_woa(:) = p_min; 
		ub_woa(:) = p_max;
		fobj_bwoa = @FBWOA_ODSTCA;
	case 'ARJOA'
		fobj_woa = @FWOA;
		lb_woa(:) = p_min;
		ub_woa(:) = p_max; 
		fobj_bwoa = @FBWOA_ARJOA;
	case 'IOJOA'
		fobj_woa = @FWOA; 
		lb_woa(:) = p_min;
		ub_woa(:) = p_max; 
		fobj_bwoa = @FBWOA_IOJOA;
	case 'OFDMA'
		fobj_woa = @FWOA; 
		lb_woa(:) = p_min;
		ub_woa(:) = p_max; 
		fobj_bwoa = @FBWOA_OFDMA;
    end

function o = FWOA(P, X) 
	mu2 = 1e20; 
	rs = 0; 
	% C5: (!!)
    fc5 = sum(P.*X, 2) - PMax; 
    flag_fc5 = fc5 > 0; 
    pnal_fc5 = sum(mu*flag_fc5.*(fc5.^2)); % === mu
    o = pnal_fc5;
	for j = 1:noSubcs
		idx = X(:, j) > 0; 
		Pj = P(:, j); 
		Hj = hArray(:, j); % ==
		Pj_off = Pj(idx); 
		Hj_off = Hj(idx); 
		no_off_users = size(Pj_off, 1); 
		if no_off_users == 0
			continue
		end 

		% interfer_j: intra-interference of users offloading to sub j
		interfer_j = zeros(no_off_users, 1); 
		for i = 1:no_off_users
			flag_less = Hj_off <= Hj_off(i); % devised in 2020/08/13
			flag_less(i) = 0; 

			% calculate the interference that each user sufferd 
			interfer_j(i) = sum(flag_less.*Pj_off.*Hj_off); 
		end

        etaj = eta(idx);
        gammaj = gamma(idx);
        pout_j = Pj_off.*Hj_off;
        % OBJECTIVE FUNCTION
        func_obj = sum((etaj + gammaj.*Pj_off)./(log2(1 + pout_j./(n0 + interfer_j))));
        % CONSTRAINT DEALING 
        % by Penalty Method 
        % C3: satisfy with lower bound condition 
        % C8: 
		
		% devised in 2020/08/12        
        fc8 = P_tol*interfer_j - pout_j; 

        % fc8 = P_tol*(interfer_j - n0) - pout_j; 
        flag_fc8 = fc8 > 0; 

        pnal_fc8 = mu2*flag_fc8.*fc8; 
        rs = rs + func_obj + sum(pnal_fc8); 
    end
   
	o = o + rs; 
end

function o = FBWOA_ODSTCA(X)
	% objective function 
	tmp = sqrt(lambda(:, 1).*f_l); 
	A = sum(X, 2); 
	fobj1 = sum(sum(lambda, 2).*A); 
	fobj3 = sum(A.*tmp)^2/f0; 
	% constraint dealing 
	fc3 = A - 1; 
	flag_fc3 = fc3 > 0; 
	pnal_fc3 = sum(mu.*flag_fc3.*(fc3.^2)); 

	o = fobj1 - fobj3 - pnal_fc3; 
end 


function o = FBWOA_ARJOA(X)
	% objective function 
	tmp = sqrt(lambda(:, 1).*f_l); 
	A = sum(X, 2); 
	fobj1 = sum(sum(lambda, 2).*A); 
	fobj3 = sum(A.*tmp)^2/f0; 
	% constraint dealing 
	fc3 = A - 1; 
	flag_fc3 = (fc3 ~= 0); 
	pnal_fc3 = sum(mu.*flag_fc3.*(fc3.^2)); 

	o = fobj1 - fobj3 - pnal_fc3; 
end 

function o = FBWOA_IOJOA(X) % A was determined  
	% objective function 
	tmp = sqrt(lambda(:, 1).*f_l); 
	A = sum(X, 2); 
	fobj1 = sum(sum(lambda, 2).*A); 
	fobj3 = sum(A.*tmp)^2/f0; 
	% constraint dealing 
	fc3 = A - Adet; 
	flag_fc3 = (fc3 ~= 0); 
	pnal_fc3 = sum(mu.*flag_fc3.*(fc3.^2)); 

	o = fobj1 - fobj3 - pnal_fc3; 
end 

function o = FBWOA_OFDMA(X)
	% objective function 
	tmp = sqrt(lambda(:, 1).*f_l); 
	A = sum(X, 2); 
	fobj1 = sum(sum(lambda, 2).*A); 
	fobj3 = sum(A.*tmp)^2/f0; 
	% constraint dealing 
	fc3 = A - 1; 
	flag_fc3 = (fc3 > 0); 
	pnal_fc3 = sum(mu.*flag_fc3.*(fc3.^2)); 

	fc4 = sum(X, 1) - 1; 
	flag_fc4 = (fc4 > 0); 
	pnal_fc4 = sum(mu.*flag_fc4.*(fc4.^2)); 

	o = fobj1 - fobj3 - pnal_fc3 - pnal_fc4; 
end 

end


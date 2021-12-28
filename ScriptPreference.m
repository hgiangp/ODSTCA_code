%-----------------------------------
% change user's preference in time
%-----------------------------------
tic
clear all
load('parameters.mat')
% noSearchAgents = 30;
% maxIter = 150;
% noSubcs = 5;
% cellRadiusMax = 250;
% cellRadiusMin = 5;
% logNormalMean = 0;
% logNormalDeviation = 8.0;

% noRealizations = 300;
% users_no = 18;
Lambda_t = 0.1:0.1:0.9;

% n0 = db2lin(-114 - 30);
% W = 1e6;

% p_min = 1e-8;
% p_max = 0.25;
% f0 = 10*1e9;
% alpha = 1*420e3;
% beta = 1000e6;

% kappa = 5e-27;
% zeta = 1;

% doTol = 1;

% mu = 1e14;
% P_tol = 1.001; % xn0 in getFunctionDetails

% f_local = 1e9*[0.5 0.8 1];
% f_user = zeros(1000, 1);
% for i = 1:1000
%     f_user(i) = f_local(randi(length(f_local), 1));
% end

po_ODSTCA = zeros(length(Lambda_t), noRealizations);
su_ODSTCA = zeros(length(Lambda_t), noRealizations);
ti_ODSTCA = zeros(length(Lambda_t), noRealizations);
en_ODSTCA = zeros(length(Lambda_t), noRealizations);

xt_lambda = cell(1, length(Lambda_t));


channelGain = zeros(users_no, noSubcs, noRealizations);
for iReal = 1:noRealizations
    [hA, ~, ~] = channelModel(users_no, noSubcs, cellRadiusMax, cellRadiusMin, logNormalMean, logNormalDeviation);
    channelGain(:, :, iReal) = hA;
end

for iK = 1:length(Lambda_t)
    lambda_t = Lambda_t(iK);
    lambda_e = 1 - lambda_t;
    lambda = [lambda_t lambda_e];
    xt_lambda(iK) = {num2str(lambda_t)};
    
    for iReal = 1:noRealizations
        [iK iReal]
        hArray = channelGain(:, :, iReal);
        t = randi(800, 1)
        f_l = f_user(t:t-1+ users_no);
        T_l = beta./f_l;
        E_l = kappa.*beta.*(f_l).^2;
        
        eta = lambda_t.*alpha./(W.*T_l);
        gamma = lambda_e.*alpha./(zeta*W.*E_l);
        
        %%%%%%%%%%%%%%%%%%%%%
        %      ODSTCA
        %%%%%%%%%%%%%%%%%%%%%
        
        T_r_ODSTCA = zeros(users_no, 1);
        E_r_ODSTCA = zeros(users_no, 1);
        Rij_ODSTCA = zeros(users_no, 1);
        
        Adetermined = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', users_no, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [~, leader_pos_bwoa, leader_pos_woa, ~, ~] = BWOA2(...
            'ODSTCA', doTol, noSearchAgents, users_no, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        % offloading vector
        A_ODSTCA = sum(leader_pos_bwoa, 2);
        off_users_no = sum(A_ODSTCA);
        
        po_ODSTCA(iK, iReal) = off_users_no/users_no;
        off_idx = (A_ODSTCA == 1);
        loc_idx = (A_ODSTCA == 0);
        Z_l_ODSTCA = loc_idx.*(lambda_t.*T_l + lambda_e.*E_l);
        
        tmp1 = sqrt(A_ODSTCA.*lambda_t.*f_l);
        f_i = f0.*(tmp1./sum(tmp1));
        
        for i = 1:users_no
            if A_ODSTCA(i) == 0
                continue;
            end
            % j: offloading subcs
            jidx = find(leader_pos_bwoa(i, :) == 1);
            % k: offloading users to subcs j that its gain smaller than i's gain
            xkj = (leader_pos_bwoa(:, jidx) > 0) & (hArray(:, jidx) < hArray(i, jidx));
            
            % data rate
            Rij_ODSTCA(i) = W*log2(1 + (leader_pos_woa(i,jidx)*hArray(i, jidx))/(n0 + ...
                sum(xkj.*leader_pos_woa(:, jidx).*hArray(:, jidx))));
            
            T_r_ODSTCA(i) = alpha/Rij_ODSTCA(i) + beta/f_i(i);
            E_r_ODSTCA(i) = leader_pos_woa(i, jidx)*alpha/(zeta*Rij_ODSTCA(i));
        end
        Z_r_ODSTCA = lambda_t.*T_r_ODSTCA + lambda_e*E_r_ODSTCA;
        
        T_i = loc_idx.*T_l + off_idx.*T_r_ODSTCA;
        E_i = loc_idx.*E_l + off_idx.*E_r_ODSTCA;
        
        Z_i = lambda_t.*(T_l - T_i)./T_l + lambda_e.*(E_l - E_i)./E_l;
        
        su_ODSTCA(iK, iReal) = sum(Z_i);
        ti_ODSTCA(iK, iReal) = sum(T_i);
        en_ODSTCA(iK, iReal) = sum(E_i);
    end
end

po_ODSTCA_1 = mean(po_ODSTCA, 2);
su_ODSTCA_1 = mean(su_ODSTCA, 2);
ti_ODSTCA_1 = mean(ti_ODSTCA, 2);
en_ODSTCA_1 = mean(en_ODSTCA, 2);

xt_lambda = cell(1, length(Lambda_t));
for i = 1:2:length(Lambda_t)
    xt_lambda(i) = {num2str(Lambda_t(i))};
end

save('preference.mat', 'Lambda_t' ,'xt_lambda','po_ODSTCA_1', 'su_ODSTCA_1', 'ti_ODSTCA_1', 'en_ODSTCA_1');

hold

ti_ODSTCA_1_aver = ti_ODSTCA_1/users_no;
en_ODSTCA_1_aver = en_ODSTCA_1/users_no;
M = length(Lambda_t)
h1 = figure(1)
hold on
yyaxis left;
plot(1:M, ti_ODSTCA_1_aver(1:M), 'b-o', 'linewidth', 2.0, 'markers', 13.0);
xlabel('User''s Preference in Time \lambda_t');
ylabel('Completion Time (s)');
yyaxis right;
ylabel('Energy Consumption(J)');
plot(1:M, en_ODSTCA_1_aver(1:M), 'b--o', 'linewidth', 2.0, 'markers', 13.0);
grid on
box on
set(gca, 'FontSize', 18, 'XLim', [1 M]);
xticks = 1:9;
set(gca, 'xtick', xticks);
set(gca, 'xticklabel',xt_lambda)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold

savefig(h1, 'prefer');

toc
t = toc
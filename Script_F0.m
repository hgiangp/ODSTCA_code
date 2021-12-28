% ---------------------------------
% change server's computing capacity
% ---------------------------------
tic
clear all

load('parameters.mat')
% noSearchAgents = 30;
% maxIter = 150;
% NoUsers = 18;
% noSubcs = 5;
% cellRadiusMax = 250;
% cellRadiusMin = 5;
% logNormalMean = 0;
% logNormalDeviation = 8.0;

% noRealizations = 300;

% lambda_t = 0.5;
% lambda_e = 1 - lambda_t;
% lambda = [lambda_t lambda_e];
% n0 = db2lin(-114 - 30);
% W = 1e6;

% p_min = 1e-8;
% p_max = 0.25;
F0 = 8:17;

% alpha = 1*420e3;
% beta = 1000e6;
% kappa = 5e-27;
% zeta = 1;

% mu = 1e14;
% P_tol = 1.001; % xn0 in getFunctionDetails

% f_local = 1e9*[0.5 0.8 1];
% f_user = zeros(1000, 1);
% for i = 1:1000
%     f_user(i) = f_local(randi(length(f_local), 1));
% end

% doTol = 1;

po_ALCA = zeros(length(F0), noRealizations);
su_ALCA = zeros(length(F0), noRealizations);

po_ARJOA = zeros(length(F0), noRealizations);
su_ARJOA = zeros(length(F0), noRealizations);

po_ODSTCA = zeros(length(F0), noRealizations);
su_ODSTCA = zeros(length(F0), noRealizations);

po_IOJOA = zeros(length(F0), noRealizations);
su_IOJOA = zeros(length(F0), noRealizations);

po_OFDMA = zeros(length(F0), noRealizations);
su_OFDMA = zeros(length(F0), noRealizations);

xt = cell(1, length(F0));

channelGain = zeros(NoUsers, noSubcs, noRealizations);
for iReal = 1:noRealizations
    [hA, ~, ~] = channelModel(NoUsers, noSubcs, cellRadiusMax, cellRadiusMin, logNormalMean, logNormalDeviation);
    channelGain(:, :, iReal) = hA;
end


for iN = 1:length(F0)
    f0 = F0(iN);
    xt(iN) = {num2str(f0)};
    f0 = f0*1e9; 

    for iReal = 1:noRealizations
        [iN iReal]
        t = randi(800, 1);
        f_l = f_user(t:t+NoUsers-1);
        T_l = beta./f_l;
        E_l = kappa.*beta.*(f_l).^2;
        
        %     %%%%%%%%%%%%%%%%%%%%%
        %     %       ALCA
        %     %%%%%%%%%%%%%%%%%%%%%
        po_ALCA(iN, iReal) = 0;
        su_ALCA(iN, iReal) = 0;

        
        eta = lambda_t.*alpha./(W.*T_l);
        gamma = lambda_e.*alpha./(zeta*W.*E_l);
        
        hArray = channelGain(:, :, iReal);
        
        %     %%%%%%%%%%%%%%%%%%%%%
        %     %       ARJOA
        %     %%%%%%%%%%%%%%%%%%%%%
                
        Adetermined = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ARJOA', NoUsers, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~] = BWOA2(...
            'ARJOA', doTol, noSearchAgents, NoUsers, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        %		cg_bwoa_ARJOA(iReal, :) = conver_curve_bwoa;
        
        % offloading vector
        A_ARJOA = sum(leader_pos_bwoa, 2);
        off_NoUsers = sum(A_ARJOA);
        
        po_ARJOA(iN, iReal) = off_NoUsers/NoUsers;
        
        su_ARJOA(iN, iReal) = leader_score_bwoa;
        
        %%%%%%%%%%%%%%%%%%%%%
        %      ODSTCA
        %%%%%%%%%%%%%%%%%%%%%
        
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', NoUsers, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~] = BWOA2(...
            'ODSTCA', doTol, noSearchAgents, NoUsers, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        % offloading vector
        A_ODSTCA = sum(leader_pos_bwoa, 2);
        off_NoUsers = sum(A_ODSTCA);
        
        po_ODSTCA(iN, iReal) = off_NoUsers/NoUsers;
        su_ODSTCA(iN, iReal) = leader_score_bwoa;        
        %%%%%%%%%%%%%%%%%%%%%
        %      IOJOA
        %%%%%%%%%%%%%%%%%%%%%
        
        Adetermined = zeros(NoUsers, 1);
        % all users independently make offloading decision
        for i = 1:NoUsers
            if rand() > 0.5
                continue
            end
            Adetermined(i) = 1;
        end
        
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('IOJOA', NoUsers, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~] = BWOA2(...
            'IOJOA', doTol, noSearchAgents, NoUsers, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        %%% ==================================
        
        % offloading vector
        
        A_IOJOA = sum(leader_pos_bwoa, 2);
        off_NoUsers = sum(A_IOJOA);
        
        po_IOJOA(iN, iReal) = off_NoUsers/NoUsers;
        su_IOJOA(iN, iReal) = leader_score_bwoa;        
        
        %%%%%%%%%%%%%%%%%%%%%
        %      OFDMA
        %%%%%%%%%%%%%%%%%%%%%
        
        Adetermined = 0; 

        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('OFDMA', NoUsers, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, conver_curve_bwoa, conver_curve_woa] = BWOA2(...
            'OFDMA', doTol, noSearchAgents, NoUsers, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        % offloading vector
        A_OFDMA = sum(leader_pos_bwoa, 2);
        off_NoUsers = sum(A_OFDMA);
        
        po_OFDMA(iN, iReal) = off_NoUsers/NoUsers;
        su_OFDMA(iN, iReal) = leader_score_bwoa;                
    end
end
po_ODSTCA = mean(po_ODSTCA, 2);
su_ODSTCA = mean(su_ODSTCA, 2);

po_ARJOA = mean(po_ARJOA, 2);
su_ARJOA = mean(su_ARJOA, 2);

po_IOJOA = mean(po_IOJOA, 2);
su_IOJOA = mean(su_IOJOA, 2);

po_OFDMA = mean(po_OFDMA, 2);
su_OFDMA = mean(su_OFDMA, 2);

po_ALCA = mean(po_ALCA,2);
su_ALCA = mean(su_ALCA,2);


save('f0com.mat', 'F0', 'xt', 'po_ALCA', 'po_ODSTCA', 'po_ARJOA', 'po_IOJOA', 'po_OFDMA', ...
    'su_ALCA', 'su_ODSTCA', 'su_ARJOA','su_OFDMA', 'su_IOJOA');

F0 = 8:17; 
h2 = figure(2)
hold on
plot(1:length(po_ODSTCA), po_ODSTCA(1:length(po_ODSTCA)), 'b-o', 'linewidth', 2.0, 'markers', 13.0);
plot(1:length(po_ARJOA), po_ARJOA(1:length(po_ARJOA)), 'g-v', 'linewidth', 2.0, 'markers', 13.0);
plot(1:length(po_IOJOA), po_IOJOA(1:length(po_IOJOA)), 'k-x', 'linewidth', 2.0, 'markers', 13.0);
plot(1:length(po_OFDMA), po_OFDMA(1:length(po_OFDMA)), 'r-s', 'linewidth', 2.0, 'markers', 13.0);
%plot(1:length(po_ALCA), po_ALCA(1:length(po_ALCA)), 'k-d', 'linewidth', 2.0, 'markers', 13.0);
grid on
box on 
set(gca, 'FontSize', 20, 'XLim', [1 length(F0)]);
xticks = 1:length(F0);
set(gca, 'xtick', xticks);
set(gca, 'xticklabel',xt)
ylabel('Offloading Percentage');
xlabel('Computing Resource f_0 (GHz)'); 
% title('N = 12', 'Interpreter', 'latex'); 
lgnd = legend({'ODSTCA', 'ARJOA', 'IOJOA', 'OFDMA'}, 'Location', 'Best'); 
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 14); 

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold

savefig(h2, 'f0com_po')

h3 = figure(3)
hold on
plot(1:length(su_ODSTCA), su_ODSTCA(1:length(su_ODSTCA)), 'b-o', 'linewidth', 2.0, 'markers', 13.0);
plot(1:length(su_ARJOA), su_ARJOA(1:length(su_ARJOA)), 'g-v', 'linewidth', 2.0, 'markers', 13.0);
plot(1:length(su_IOJOA), su_IOJOA(1:length(su_IOJOA)), 'k-x', 'linewidth', 2.0, 'markers', 13.0);
plot(1:length(su_OFDMA), su_OFDMA(1:length(su_OFDMA)), 'r-s', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(su_ALCA), su_ALCA(1:length(su_ALCA)), 'k-d', 'linewidth', 2.0, 'markers', 13.0);
grid on
box on 
set(gca, 'FontSize', 17.5, 'XLim', [1 length(F0)]);
xticks = 1:length(F0);
set(gca, 'xtick', xticks);
set(gca, 'xticklabel',xt)
ylabel('System Utility');
xlabel('Computing Resource f_0 (GHz)');
lgnd = legend({'ODSTCA', 'ARJOA', 'IOJOA', 'OFDMA'}, 'Location', 'Best'); 
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 17.5); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold
savefig(h3, 'f0com_su')

toc


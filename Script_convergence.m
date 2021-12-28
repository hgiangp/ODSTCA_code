%-----------------------------------
% check number of iterations until the algorithm is converged
%-----------------------------------
tic
clear all
load('parameters.mat')
% noSearchAgents = 30;
maxIter = 500;
NoUsers = [5 15 25];
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

% p_min = 1e-6;
% p_max = 0.25;
% f0 = 10*1e9;
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

doTol = 0;

cv_bwoa = zeros(noRealizations, maxIter, length(NoUsers)); 
cv_woa = zeros(noRealizations, maxIter, length(NoUsers)); 

cv_bwoa_mean = zeros(length(NoUsers), maxIter); 
cv_woa_mean = zeros(length(NoUsers), maxIter); 

xt = cell(1, length(NoUsers));

for iN = 1:length(NoUsers)
    users_no = NoUsers(iN);
    xt(iN) = {num2str(users_no)};
    channelGain = zeros(users_no, noSubcs, noRealizations);
    for iReal = 1:noRealizations
        [hA, ~, ~] = channelModel(users_no, noSubcs, cellRadiusMax, cellRadiusMin, logNormalMean, logNormalDeviation);
        channelGain(:, :, iReal) = hA;
    end
    
    for iReal = 1:noRealizations
        [iN iReal]
        t = randi(800, 1);
        f_l = f_user(t:t+users_no-1);
        T_l = beta./f_l;
        E_l = kappa.*beta.*(f_l).^2;
                
        eta = lambda_t.*alpha./(W.*T_l);
        gamma = lambda_e.*alpha./(zeta*W.*E_l);
        
        hArray = channelGain(:, :, iReal);
                
        %%%%%%%%%%%%%%%%%%%%%
        %      ODSTCA
        %%%%%%%%%%%%%%%%%%%%%
        
        Adetermined = 1; 
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', users_no, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, conver_curve, conver_curve_woa] = BWOA2(...
            'ODSTCA', doTol, noSearchAgents, users_no, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        cv_bwoa(iReal, :, iN) = conver_curve; 
        cv_woa(iReal, :, iN) = conver_curve_woa; 
    end
    cv_bwoa_mean(iN, :) = mean(cv_bwoa(:, :, iN)); 
    cv_woa_mean(iN, :) = mean(cv_woa(:, :, iN)); 
end

save('conver.mat','NoUsers', 'xt', 'cv_bwoa', 'cv_woa', 'cv_bwoa_mean', 'cv_woa_mean')

maxIter = 150
xt = {'10', '25', '50', '100', '150'}; 

h1 = figure(1)
hold on
plot(1:maxIter, cv_bwoa_mean(1, 1:maxIter), 'b-', 'linewidth', 2.0, 'markers', 13.0);
plot(1:maxIter, cv_bwoa_mean(2, 1:maxIter), 'b--', 'linewidth', 2.0, 'markers', 13.0);
plot(1:maxIter, cv_bwoa_mean(3, 1:maxIter), 'b:', 'linewidth', 2.0, 'markers', 13.0);
grid on 

set(gca, 'FontSize', 17.5, 'XLim', [1 maxIter]);
xticks = [10, 25, 50, 100, 150]; 
set(gca, 'xtick', xticks);
set(gca, 'xticklabel',xt); 

xlabel('Iteration Index');
ylabel('System Utility');
lgnd = legend({'ODSTCA, N = 5', 'ODSTCA, N = 15', 'ODSTCA, N = 25'}, 'Location', 'Best')
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 17.5); 
box on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold

savefig(h1, 'conver'); 

toc

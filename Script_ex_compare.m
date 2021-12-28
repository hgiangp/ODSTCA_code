%------------------------------
% compare ODSTCA with exhaustive search in term of convergence behavior and runtime
% ----------------------------- 
load('parameters.mat')
% noSearchAgents = 30;
% maxIter = 150;
NoUsers = 2:7;
% noSubcs = 3;
% cellRadiusMax = 250;
% cellRadiusMin = 5;
% logNormalMean = 0;
% logNormalDeviation = 8.0;

% noRealizations = 1;

% lambda_t = 0.5;
% lambda_e = 1 - lambda_t;
% lambda = [lambda_t lambda_e];
% n0 = db2lin(-114 - 30);
% W = 1e6;

% p_min = 1e-8;
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

% doTol = 1;

po_ODSTCA = zeros(length(NoUsers), noRealizations);
su_ODSTCA = zeros(length(NoUsers), noRealizations);
time_ODSTCA = zeros(length(NoUsers), noRealizations);

po_EX = zeros(length(NoUsers), noRealizations);
su_EX = zeros(length(NoUsers), noRealizations);
time_EX = zeros(length(NoUsers), noRealizations);

for iN = 1:length(NoUsers)
    users_no = NoUsers(iN);
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
        
        Adet = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', users_no, noSubcs, hArray, p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adet);
        
        % Exhaustive search
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, time] = exhaustive(noSearchAgents, users_no, noSubcs, maxIter, lb_woa, ub_woa, fobj_woa, lambda, f_l, f0);
        po_EX(iN, iReal) = sum(sum(leader_pos_bwoa))/users_no;
        su_EX(iN, iReal) = leader_score_bwoa;
        time_EX(iN, iReal) = time;
        % ODSTCA
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~, ~, time] = BWOA2(...
            'ODSTCA', doTol, noSearchAgents, users_no, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adet);
        po_ODSTCA(iN, iReal) = sum(sum(leader_pos_bwoa))/users_no;
        su_ODSTCA(iN, iReal) = leader_score_bwoa;
        time_ODSTCA(iN, iReal) = time;
    end
end
po_ODSTCA = mean(po_ODSTCA, 2);
su_ODSTCA = mean(su_ODSTCA, 2);
time_ODSTCA = mean(time_ODSTCA, 2);

po_EX = mean(po_EX, 2);
su_EX = mean(su_EX, 2);
time_EX = mean(time_EX, 2);
NoUsers = 2:8;
xt = {'2' ,'3', '4', '5', '6', '7', '8'};
save('exhastive.mat', 'NoUsers', 'xt', 'po_ODSTCA', 'su_ODSTCA', 'time_ODSTCA', 'po_EX', 'su_EX', 'time_EX');
xticks = 2:8;


h2 = figure(2)
hold on
plot(NoUsers, su_ODSTCA, 'b-o', 'linewidth', 2.0, 'markers', 13.0);
plot(NoUsers, su_EX, 'r-s', 'linewidth', 2.0, 'markers', 13.0);
grid on
set(gca, 'FontSize', 17.5, 'XLim', [2 7]);
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xt);
xlabel('Number of Users');
ylabel('System Utility');
legend({'ODSTCA', 'EX'}, 'Location', 'Best', 'FontSize', 17.5);
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
x = [0.3 0.5], y = [0.6 0.5]
dim = [0.2 0.2 0.1 0.1]
annotation('ellipse',dim)
a = annotation('textarrow',x,y,'String','Optimality gap = 1.22%', 'FontSize', 24);
savefig(h2, 'ex_su');

h3 = figure(3)
ax(1) = axes();
l1=line('parent',ax(1),'xdata',NoUsers,'ydata',time_ODSTCA, 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 13.0, 'Color', 'b');
set(ax(1), 'XLim', [2 7], 'FontSize', 17.5);
set(ax(1), 'xtick', xticks);
set(ax(1), 'xticklabel',xt);
xlabel('parent', ax(1), 'Number of Users');
ylabel('parent', ax(1), 'Algorithm Runtime (s)');
grid on
box on

ax = ax(1);
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

ax(2) = axes('position', [0.18 0.57 0.25 0.4]);
l2=line('parent',ax(2),'xdata',NoUsers,'ydata',time_EX,'LineWidth', 2.0, 'Marker', 's', 'MarkerSize', 13.0, 'Color', 'r');
set(ax(2), 'XLim', [2 7], 'xticklabel', [], 'FontSize', 17.5);
axis tight
box on
grid on

lngd = legend( [l1;l2] , {'ODSTCA','EX'}, 'FontSize', 20);
savefig(h3, 'ex_time');


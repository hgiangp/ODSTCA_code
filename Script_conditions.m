%------------------------------------
% compare ODSTCA in 3 version when changing the number of users: N
% 1. Not check any condition
% 2. Check only condition 1
% 3. Check both conditions
%------------------------------------
tic
clear all

load('parameters.mat')

% noSearchAgents = 30;
% maxIter = 150;
NoUsers = [4 5 6:2:24];
% noSubcs = 5;
% cellRadiusMax = 250;
% cellRadiusMin = 5;
% logNormalMean = 0;
% logNormalDeviation = 8.0;
% noRealizations = 150;
% lambda_t = 0.5;
% lambda_e = 1 - lambda_t;
% lambda = [lambda_t lambda_e];
% n0 = db2lin(-100 - 30);
% W = 1e6;

% p_min = 1e-6;
% p_max = 0.25;
% f0 = 10*1e9;
% alpha = 1*420e3;
% beta = 2000e6;
% kappa = 5e-27;
% zeta = 1;

% mu = 1e14;
% P_tol = 0.01; % xn0 in getFunctionDetails

% f_local = 1e9*[0.5 0.8 1];
% f_user = zeros(1000, 1);
% for i = 1:1000
%     f_user(i) = f_local(randi(length(f_local), 1));
% end

% doTol = 1;

su_ODSTCA_ot = zeros(length(NoUsers), noRealizations); %optimized time
no_WOA_ot = zeros(length(NoUsers), noRealizations);
time_ot = zeros(length(NoUsers), noRealizations);

su_ODSTCA_nm = zeros(length(NoUsers), noRealizations); % normal
no_WOA_nm = zeros(length(NoUsers), noRealizations);
time_nm = zeros(length(NoUsers), noRealizations);

su_ODSTCA_first = zeros(length(NoUsers), noRealizations); %first version, without checking any conditions
no_WOA_first = zeros(length(NoUsers), noRealizations);
time_first = zeros(length(NoUsers), noRealizations);

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
        %   Condition 2
        %%%%%%%%%%%%%%%%%%%%%
        
        Adetermined = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', users_no, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~, no_WOA_run, time] = BWOA2(...
            'ODSTCA', doTol, noSearchAgents, users_no, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
            gamma, eta, hArray, n0, p_min, p_max, Adetermined);
        
        su_ODSTCA_ot(iN, iReal) = leader_score_bwoa;
        
        no_WOA_ot(iN, iReal) = no_WOA_run;
        time_ot(iN, iReal) = time;
        
        
        %%%%%%%%%%%%%%%%%%%%%
        %   Condition 1
        %%%%%%%%%%%%%%%%%%%%%
        
        
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', users_no, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~, no_WOA_run, time] = BWOA1(...
            'ODSTCA', doTol, noSearchAgents, users_no, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa);
        
        su_ODSTCA_nm(iN, iReal) = leader_score_bwoa;
        
        no_WOA_nm(iN, iReal) = no_WOA_run;
        time_nm(iN, iReal) = time;
        
        %%%%%%%%%%%%%%%%%%%%%
        %   No condition
        %%%%%%%%%%%%%%%%%%%%%
        Adetermined = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails('ODSTCA', users_no, noSubcs, hArray, ...
            p_min, p_max, P_tol, n0, eta, gamma, mu, lambda, f0, f_l, Adetermined);
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~, no_WOA_run, time] = BWOA(...
            'ODSTCA', doTol, noSearchAgents, users_no, noSubcs, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa);
        
        su_ODSTCA_first(iN, iReal) = leader_score_bwoa;
        
        no_WOA_first(iN, iReal) = no_WOA_run;
        time_first(iN, iReal) = time;
    end
end

su_ODSTCA_ot = mean(su_ODSTCA_ot, 2);
no_WOA_ot = mean(no_WOA_ot, 2);
time_ot = mean(time_ot, 2);

su_ODSTCA_nm = mean(su_ODSTCA_nm, 2);
no_WOA_nm = mean(no_WOA_nm, 2);
time_nm = mean(time_nm, 2);

su_ODSTCA_first = mean(su_ODSTCA_first, 2);
no_WOA_first = mean(no_WOA_first, 2);
time_first = mean(time_first, 2);

save('script_time.mat', 'su_ODSTCA_ot', 'no_WOA_ot', 'time_ot', 'su_ODSTCA_nm', 'no_WOA_nm', 'time_nm', 'su_ODSTCA_first', 'no_WOA_first', 'time_first');

NoUsers = [4 5 6:2:24];

xticks = [4, 5, 8:4:24]; 
xt = cell(1, length(xticks)); 
for i = 1:length(xticks) 
    xt(i) = {num2str(xticks(i))}; 
end


h1 = figure(1)
hold on
plot(NoUsers, su_ODSTCA_ot(1:length(su_ODSTCA_ot)), 'b-o', 'linewidth', 2.0, 'markers', 13.0);
plot(NoUsers, su_ODSTCA_nm(1:length(su_ODSTCA_nm)), 'g-v', 'linewidth', 2.0, 'markers', 13.0);
plot(NoUsers, su_ODSTCA_first(1:length(su_ODSTCA_first)), 'r-s', 'linewidth', 2.0, 'markers', 13.0);
grid on
set(gca, 'XLim', [4 24], 'FontSize', 20); 
set(gca, 'xtick', xticks);
set(gca, 'xticklabel',xt);
xlabel('Number of Users');
ylabel('System Utility');
lgnd = legend({'ODSTCA', 'Con. 1', 'No Con'}, 'FontSize', 14, 'Location', 'Best')
box on
hold 
savefig(h1, 'time_su.fig')

h2 = figure(2)
ax(1) = axes();  
l1=line('parent',ax(1),'xdata',NoUsers,'ydata',no_WOA_ot, 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 13.0, 'Color', 'b');
l2=line('parent',ax(1),'xdata',NoUsers,'ydata',no_WOA_nm, 'LineWidth', 2.0, 'Marker', 'v', 'MarkerSize', 13.0, 'Color', 'g');
set(ax(1), 'XLim', [4 24], 'FontSize', 20); 
set(ax(1), 'xtick', xticks);
set(ax(1), 'xticklabel',xt);
xlabel('parent', ax(1), 'Number of Users');
ylabel('parent', ax(1), 'Number of WOA Calls');
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

ax(2) = axes('position', [0.7 0.30 0.25 0.4]);
l3=line('parent',ax(2),'xdata',NoUsers,'ydata',no_WOA_first,'LineWidth', 2.0, 'Marker', 's', 'MarkerSize', 13.0, 'Color', 'r');
set(ax(2), 'XLim', [4 24], 'xticklabel', [], 'FontSize', 14.0)
box on 
grid on 
axis tight 
lngd = legend([l1,l2, l3], {'ODSTCA','Con. 1','No. Con'}, 'FontSize', 14);
savefig(h2, 'time_WOA.fig')


h3 = figure(3)
ax(1) = axes();  
l1=line('parent',ax(1),'xdata',NoUsers,'ydata',time_ot, 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 13.0, 'Color', 'b');
l2=line('parent',ax(1),'xdata',NoUsers,'ydata',time_nm, 'LineWidth', 2.0, 'Marker', 'v', 'MarkerSize', 13.0, 'Color', 'g');
set(ax(1), 'XLim', [4 24], 'FontSize', 20); 
set(ax(1), 'xtick', xticks);
set(ax(1), 'xticklabel',xt);
xlabel('parent', ax(1), 'Number of Users');
ylabel('parent', ax(1), 'BWOA Run Time (s)');
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
l3=line('parent',ax(2),'xdata',NoUsers,'ydata',time_first,'LineWidth', 2.0, 'Marker', 's', 'MarkerSize', 13.0, 'Color', 'r');
set(ax(2), 'XLim', [4 24], 'xticklabel', [], 'FontSize', 14.0)
axis tight 
box on 
grid on 

lngd = legend( [l1;l2;l3] , {'ODSTCA','Con. 1','No. Con'}, 'FontSize', 14);
savefig(h3, 'time_ti.fig')
toc

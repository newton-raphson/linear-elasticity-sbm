clear;clc;

ld = 3; % linewidth
ft = 30;

%L2: 6.33878e-5(5)->9.77857e-6(7)

X=[1/2^5,1/2^6,1/2^7,1/2^8,1/2^9];

Y_wo=[1.02476e-5,4.3832e-6,1.85174e-6,9.7413e-7,5.19755e-7]; % 問題

% L2 - domain
Y_w=[0.00161458,0.000398943,0.000104486,2.66728e-5,6.97317e-6];
Y_nit=[0.0017201,0.000655008,0.000196168,6.53007e-5,4.11446e-5]; % small code


% domain
Y_w3=[0.00249289,0.000474527,0.000111275,2.77274e-05,7.15476e-6];
Y_nit3=[0.00367813,0.00160729,0.000758368,0.000264413,0.000170352]; % small


%%
% Y_w=Y_w3;
% Y_nit=Y_nit3;

%%

slope_W=cal_slope(X,Y_w);
slope_NIT=cal_slope(X,Y_nit);


%%
figure('color','w','position',[300,300,1000,600]); % 寬/高

loglog(X,Y_w,'-s','LineWidth',ld,'color','r');hold on;
loglog(X,Y_nit,'-s','LineWidth',ld,'color','b');hold on;

xlabel('h','FontSize',ft,'FontName','Times')
ylabel('||u-u_{h}||','FontSize',ft,'FontName','Times')
set(gca,'FontSize',ft,'FontName','Times New Roman');
h=legend({['L2 norm (SBM):',num2str(slope_W)],['L2 norm (nitsche):',num2str(slope_NIT)]},'Location','northeastoutside');

set(h,'box','off');
grid on
grid minor

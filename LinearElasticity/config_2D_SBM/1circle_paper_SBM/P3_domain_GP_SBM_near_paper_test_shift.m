clear;clc;

co = [
0 0 0; 
0.137254902   0.254901961	 1.0; %blue
0 0.6 0.15; %green  
.5 0 .5;
0.350000 0.000000 1.000000;
0.760000 0.700000 0.500000;
0 0.5 0;
1 0 0;
0 0.447 0.741
];


set(groot,'defaultAxesColorOrder',co);

ld = 3; % linewidth
ft = 30;

%L2: 6.33878e-5(5)->9.77857e-6(7)

X=[1/2^5,1/2^6,1/2^7,1/2^8,1/2^9];

% L2 - domain
Y_w1=[0.00142931,0.000379867,0.000100732,2.63925e-05,7.08174e-06];
Y_nit1=[0.00151807,0.000556843,0.000185507,9.06167e-5,2.90977e-5]; % small code

%2 domain (shortest)
Y_w1_short=[0.00140632,0.000366972,9.65028e-05,2.52762e-05 ,6.61649e-06];



% L2 - domain
Y_w=[0.00161458,0.000398943,0.000104486,2.66728e-5,6.97317e-6];
Y_nit=[0.0017201,0.000655008,0.000196168,6.53007e-5,4.11446e-5]; % small code

%2 domain (shortest)
Y_w_short=[0.00187384 ,0.000419524 ,0.000103149,2.61453e-05 ,6.7549e-06];

% domain
Y_w3=[0.00249289,0.000474527,0.000111275,2.77274e-05,7.15476e-6];
Y_nit3=[0.00367813,0.00160729,0.000758368,0.000264413,0.000170352]; % small

%3 domain (shortest)
Y_w3_short=[0.00358614,0.000572508,0.000127545,2.91586e-05,7.53695e-06];

Y_w3_short_shift=[0.00337927,0.000776463,0.000228003,4.95418e-5,1.3355e-5];

Y_w1_short_shift=[0.00210489,0.000560686,0.000148689,3.87123e-5,1.011881e-5];

slope_short1=cal_slope(X,Y_w1_short);
slope_short1_shift=cal_slope(X,Y_w1_short_shift);

%%
% figure('color','w','position',[300,300,1000,600]); % 寬/高
figure('color','w','position',[300,300,1400,800]); % 寬/高

loglog(X,Y_w1_short,'-s','LineWidth',ld);hold on;
loglog(X,Y_w1_short_shift,'-s','LineWidth',ld);hold on;

xlabel('h','FontSize',ft,'FontName','Times')
ylabel('||u-u_{h}||','FontSize',ft,'FontName','Times')
set(gca,'FontSize',ft,'FontName','Times New Roman');
h=legend( ...
    {    ['L2 norm:',num2str(slope_short1)] , ...
     ['L2 norm (shift the center and expand the domain):',num2str(slope_short1_shift)] ...
    },'Location','northeastoutside');

set(h,'box','off');
grid on
grid minor

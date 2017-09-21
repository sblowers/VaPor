Tt(~GM_WM) = NaN;

% delta = 0.0031;
delta = 0.0031;
Tm = 0.19;
R = 0.012;
brain_vector = fliplr(linspace(VoxelSize,R,40));
TSE_function = BloodTemp + Tm - (H_Out*delta*(BloodTemp+Tm-T_Out)/(mean(Kc(GM_WM))*besseli(1,R/delta)+H_Out*delta*besseli(0,R/delta))).*besseli(0,(R-brain_vector)/delta);


% VaPor_function = squeeze(Ts(30,36,squeeze(~isnan(Ts(30,36,:)))));
% VaPor_function = VaPor_function(end-39:end);
% VaPor_function2 = squeeze(Ts(30,30,squeeze(~isnan(Ts(30,30,:)))));
% VaPor_function2 = VaPor_function2(end-39:end);
% VaPor_function3 = squeeze(Ts(30,24,squeeze(~isnan(Ts(30,24,:)))));
% VaPor_function3 = VaPor_function3(end-39:end);

VaPor_function = squeeze(Tt(30,36,squeeze(~isnan(Tt(30,36,:)))));
VaPor_function = VaPor_function(end-39:end);
VaPor_function2 = squeeze(Tt(30,30,squeeze(~isnan(Tt(30,30,:)))));
VaPor_function2 = VaPor_function2(end-39:end);
VaPor_function3 = squeeze(Tt(30,24,squeeze(~isnan(Tt(30,24,:)))));
VaPor_function3 = VaPor_function3(end-39:end);

figure
hold on
plot(linspace(VoxelSize,40*VoxelSize,40),VaPor_function,'k-o','MarkerFaceColor','black')
plot(linspace(VoxelSize,40*VoxelSize,40),VaPor_function2,'k-s','MarkerFaceColor','black')
plot(linspace(VoxelSize,40*VoxelSize,40),VaPor_function3,'k-^','MarkerFaceColor','black')
plot(linspace(VoxelSize,40*VoxelSize,40),TSE_function,'k--o')

grid on

xlabel('Radius [m]')
ylabel('Temperature [^{o}C]')
legend('VaPor Loc-1','VaPor Loc-2','VaPor Loc-3','TSE','Location','southwest')
% legend('Pennes Loc-1','Pennes Loc-2','Pennes Loc-3','TSE','Location','southwest')

set(gca,'FontSize',15)

% % VaPor_function_met2 = squeeze(Ts(30,36,squeeze(~isnan(Ts(30,36,:)))));

% figure
% hold on
% plot(linspace(h,40*h,40),CoVat_function,'k-o','MarkerFaceColor','black')
% plot(linspace(h,40*h,40),CoVat_function_met1,'k-s','MarkerFaceColor','black')
% plot(linspace(h,40*h,40),CoVat_function_met2,'k-^','MarkerFaceColor','black')
% plot(linspace(h,40*h,40),TSE_function,'k--o')
% plot(linspace(h,40*h,40),TSE_function_met1,'k--s')
% plot(linspace(h,40*h,40),TSE_function_met2,'k--^')
% 
% legend('VaPor-FET T_{m}=0.19^{o}C','VaPor-FET T_{m}=0.57^{o}C','VaPor-FET T_{m}=0^{o}C','TSE T_{m}=0.19^{o}C','TSE T_{m}=0.57^{o}C','TSE T_{m}=0^{o}C','Location','southwest')
% 
% grid on
% 
% xlabel('Radius - m')
% ylabel('Temperature - ^{o}C')
% 
% set(gca,'FontSize',15)
% 
% figure
% hold on
% 
% % % VaPor_function_perf1= squeeze(Ts(30,36,squeeze(~isnan(Ts(30,36,:)))));
% 
% plot(linspace(h,40*h,40),CoVat_function,'k-o','MarkerFaceColor','black')
% plot(linspace(h,40*h,40),CoVat_function_met1,'k-s','MarkerFaceColor','black')
% plot(linspace(h,40*h,40),CoVat_function_perf2,'k-^','MarkerFaceColor','black')
% plot(linspace(h,40*h,40),TSE_function,'k--o')
% plot(linspace(h,40*h,40),TSE_function_perf1,'k--s')
% plot(linspace(h,40*h,40),TSE_function_perf2,'k--^')
% 
% legend('VaPor-FET CBF=81ml/100g/min','VaPor-FET CBF=121.5ml/100g/min','VaPor-FET CBF=40.5ml/100g/min','TSE CBF=81ml/100g/min','TSE CBF=121.5ml/100g/min','TSE CBF=40.5ml/100g/min','Location','southwest')
% 
% grid on
% 
% xlabel('Radius - m')
% ylabel('Temperature - ^{o}C')
% 
% set(gca,'FontSize',15)
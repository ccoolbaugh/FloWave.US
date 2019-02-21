function FinalPlots(AllPSVTime,AllPSV,AllPDVTime,AllPDV,AllEDVTime,AllEDV,AllMBF1)

% Plot Cardiac Cycle Points 
figure, grid on, hold on;
plot(AllPSVTime,AllPSV, '--go','DisplayName', 'Peak Systole', 'LineWidth', 1.1, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 4); hold all;
plot(AllPDVTime,AllPDV,'--bo','DisplayName','Nadir Diastole','LineWidth',1.1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize', 4);hold all;
plot(AllEDVTime,AllEDV,'--mo','DisplayName','End Diastole','LineWidth',1.1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 4);hold off;
title('Blood Flow vs. Time','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Blood Flow (ml/min)','FontSize',14)
legend('Peak Systole','Diastole','End Diastole')

% Plot Mean Blood Flow
figure, axis tight; 
plot(AllPSVTime,AllMBF1,'k','DisplayName','Mean Blood Flow','LineWidth', 1.1, 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 4), grid on;
xlabel('Time (s)', 'FontSize', 14),
ylabel('Mean Blood Flow (ml/min', 'FontSize', 14),

end


function []=PlotMCAVTerror(M_MC,MCerror,MCAVTerror)
% plot in order to compare the error computed by the Antithetic Variables
% Technique and the MonteCarlo error
%
%INPUT
% M_MC:       row vector number of simulations   
% MCerror:    row vector of MonteCarlo errors    
% MCAVTerror: row vector of errors computed with Antithetic Variables Technique 

% we want to verify that the error computed with the Antithetic Variable
% Technique is better than the Monte-Carlo one
figure()
grid on;
loglog(M_MC,MCerror,'r','LineWidth',2);
hold on 
loglog(M_MC,MCAVTerror,'b','LineWidth',2);
loglog(M_MC,1./sqrt(M_MC),'g','LineWidth',2);
xlabel('M simulations');
legend('MCerror','MCAVTerror','1/sqrt(M)');
hold off
disp(['Error improvement factor for the largest M: ' num2str(MCAVTerror(end)/MCerror(end))]);
end % function PlotMCAVTerror

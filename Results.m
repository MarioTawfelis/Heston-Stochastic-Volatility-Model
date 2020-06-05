% This script will execute Monte Carlo and Integration scripts to view results

MonteCarlo_Live  % Run Monte Carlo Live script

Integration_Live % Run Integration Live script

%%
disp(' ')
fprintf('Method                               Call        Put\n')
fprintf('------------------------------------------------------\n')
fprintf('Monte Carlo                         %5.4f      %5.4f   \n', MC_call,MC_put);
fprintf('with standard deviation of          %5.4f      %5.4f   \n', MC_call_stdev,MC_put_stdev);
disp(' ')
fprintf('Schoutens et al. (2004):\n')
fprintf('          Gauss-Laguerre            %5.4f      %5.4f \n', SchoutensCall1,SchoutensPut1);
fprintf('          Trapz                     %5.4f      %5.4f \n', SchoutensCall2,SchoutensPut2);
fprintf('Cui et al. (2017):\n')
fprintf('          Gauss-Laguerre            %5.4f      %5.4f \n', CuiCall1,CuiPut1);
fprintf('          Trapz                     %5.4f      %5.4f \n', CuiCall2,CuiPut2);
fprintf('------------------------------------------------------\n')

disp(' ')

fprintf('Absolute Error                       Call        Put\n')
fprintf('------------------------------------------------------\n')
fprintf('Schoutens-Cui:                      %1.4f      %1.4f\n', SchoutensCall1-CuiCall1,SchoutensPut1-CuiPut1);
fprintf('Schoutens-Monte Carlo               %5.4f      %5.4f\n', SchoutensCall1-MC_call,SchoutensPut1-MC_put);
fprintf('Monte Carlo-Cui:                    %1.4f      %1.4f\n', MC_call-CuiCall1,MC_put-CuiPut1);
fprintf('------------------------------------------------------\n')

disp(' ')

fprintf('Execution time\n')
fprintf('------------------------------------------------------\n')
fprintf("Monte Carlo:\n")
fprintf("\tInitialisation execution time   %3.1f seconds\n",MCInitialisation)
fprintf("\tSimulation execution time       %3.1f seconds\n",MCSimulation)
fprintf("\tTotal execution time            %3.1f seconds\n",MCInitialisation+MCSimulation)
fprintf("\nNumerical methods:\n")
fprintf("\tInitialisation execution time    %3.1f seconds\n",NMInitialisation)
fprintf("\tSimulation execution time        %3.1f seconds\n",NMSimulation)
fprintf("\tTotal execution time             %3.1f seconds\n",NMInitialisation+NMSimulation)
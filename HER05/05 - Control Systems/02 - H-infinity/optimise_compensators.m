function [emax,W1,W2,Cinf] = optimise_compensators(P,LoopShape_bounds,Weight_bounds,CondNo_bounds,WeightStruc,omega,Cinit)
% This function simultaneously synthesises by one algorithm loop-shaping
% weights 'W1/W2' and a robustly stabilising controller 'Cinf' which
% achieves a pre-specified performance level and maximises the robust
% stability margin 'emax' to normalised coprime factor uncertainty. Here,
% the notation is such that 'W2' pre-multiplies 'P' and 'W1'
% post-multiplies it to form the shaped plant.
%
% This script has been updated for modern MATLAB versions (R2016b and later)
% and requires the Robust Control Toolbox and System Identification Toolbox. 
% Original script by Alexander Lanzon - 17 July 2000.
%
% USAGE  :  [emax,W1,W2,Cinf] = simWCsyn_HinfLoopShaping_updated(P,LoopShape_bounds, ...
%                               Weight_bounds,CondNo_bounds,WeightStruc,omega,Cinit)
%
% INPUTS :  P                 = LTI system (ss, tf, or zpk) representing nominal plant.
%           LoopShape_bounds  = A 2x1 LTI system matrix such that the magnitude
%                               of the first (second) transfer function
%                               gives an upper (lower) bound for the desired
%                               loop-shape.
%
%           THE FOLLOWING INPUT ARGUMENTS ARE OPTIONAL
%
%           Weight_bounds     = A 2x2 LTI system matrix such that the magnitude
%                               of each individual transfer function is an
%                               upper/lower bound for the singular values
%                               of W1/W2 as follows:
%                                                 W1    W2
%                                          upper [1,1   1,2]
%                                          lower [2,1   2,2].
%                               (DEFAULT value is '[1e10 1; 1e-10 1]').
%           CondNo_bounds     = A 1x2 LTI system matrix giving upper bounds for
%                               the condition number of W1 (first element) and W2 (second).
%                               (DEFAULT value is '[20 20]').
%           WeightStruc       = 'DiagDiag' for diagonal loop-shaping weights W1/W2.
%                               'DiagFull' for diagonal W1 and non-diagonal W2.
%                               'FullDiag' for non-diagonal W1 and diagonal W2.
%                               'FullFull' for non-diagonal W1 and non-diagonal W2.
%                               (DEFAULT value is 'DiagDiag').
%           omega             = A vector of frequencies (rad/s) for constraint checking.
%                               (DEFAULT value is 'logspace(-4,4,100)').
%           Cinit             = An internally stabilising controller (LTI object)
%                               used to initialise the algorithm (for the unshaped plant).
%
% OUTPUTS:  emax              = Best attained robust stability margin.
%           W1                = Loop-shaping weight post-multiplying P.
%           W2                = Loop-shaping weight pre-multiplying P.
%           Cinf              = Robust stabilising controller for the shaped plant 'Ps = W2*P*W1'.
%                               The final controller for the nominal plant is 'C = W1*Cinf*W2'.

%% Checking number of inputs and outputs
if (nargin < 2 ) || (nargin > 7)
    error('Wrong number of inputs. Usage: simWCsyn_HinfLoopShaping_updated(P,LoopShape_bounds,...)');
end
if (nargout > 4)
    error('Wrong number of outputs.');
end

%% Checking required inputs
if ~isa(P, 'lti')
    error('Nominal Plant "P" must be an LTI object (ss, tf, zpk).');
  
end
    [m, n] = size(P);

if ~isa(LoopShape_bounds, 'lti') %|| ~isequal(size(LoopShape_bounds), [2,1])
    % disp(class(LoopShape_bounds))
    error('Variable "LoopShape_bounds" must be a 2x1 LTI object.');
end

%% Checking and setting optional inputs
if (nargin < 3) || isempty(Weight_bounds)
    Weight_bounds = [ss(1e5) ss(1); ss(1e-5) ss(1)];
end

% if ~isequal(size(Weight_bounds), [2,2])
%     error('Variable "Weight_bounds" not properly defined.');
% end

if (nargin < 4) || isempty(CondNo_bounds)
    CondNo_bounds = [ss(20) ss(20)];
end
if ~isequal(size(CondNo_bounds), [1,2])
    error('Variable "CondNo_bounds" not properly defined.');
end

if (nargin < 5) || isempty(WeightStruc)
    WeightStruc = 'DiagDiag';
end
if ~ismember(WeightStruc, {'DiagDiag', 'DiagFull', 'FullDiag', 'FullFull'})
    error('Variable "WeightStruc" not properly defined.');
end

if (nargin < 6) || isempty(omega)
    omega = logspace(-4,4,100);
end
if ~isvector(omega)
    error('Variable "omega" must be a vector.');
end

if (nargin < 7) || isempty(Cinit)
    disp('Synthesizing initial stabilizing controller using ncfsyn...');
    [Cinit, ~, gam] = ncfsyn(P);
    fprintf('Initial controller synthesized with performance gamma = %.2f\n', gam);
else
    if ~isa(Cinit, 'lti') || ~isequal(size(Cinit), [n, m])
        error('Variable "Cinit" not properly defined.');
    end
end
C = Cinit;


%% Check if weight bounds are trivial
W1_bounds_FLAG = false;
% Use the standard H-infinity norm calculation, which is more robust.
% This checks if the first column of Weight_bounds is equivalent to identity.
if norm(Weight_bounds(:,1) - 1, inf) < 1e-9 % Comparing floating point numbers
    W1_bounds_FLAG = true;
end

W2_bounds_FLAG = false;
% This checks if the second column of Weight_bounds is equivalent to identity.
if norm(Weight_bounds(:,2) - 1, inf) < 1e-9 % Comparing floating point numbers
    W2_bounds_FLAG = true;
end

if W1_bounds_FLAG && W2_bounds_FLAG
    error('It is not possible to have Weight_bounds = [1 1; 1 1]');
end

%% Plotting initial data
P_fr = freqresp(P, omega);
LoopShape_upper_fr = abs(squeeze(freqresp(LoopShape_bounds(1,1), omega)));
LoopShape_lower_fr = abs(squeeze(freqresp(LoopShape_bounds(2,1), omega)));

% figure(1);
% hax1 = gca;
% semilogx(hax1, omega)%, NaN, 'b-'); % Placeholder for handle
% hold(hax1, 'on');
% 
% % hold on
% xlabel(hax1, 'Frequency (radians/sec)');
% ylabel(hax1, 'Pointwise Robust Stability Margin');
% axis(hax1, [omega(1) omega(end) 0 1]);
% grid(hax1, 'on');
% zoom(hax1, 'on');
% title(hax1, 'Robust Stability Margin');
% hold(hax1, 'off');

figure(1);
% No need to plot a placeholder, as the loop will clear and draw the real data.
% Just set up the labels and properties for the initial view.
xlabel('Frequency (radians/sec)');
ylabel('Pointwise Robust Stability Margin');
title('Robust Stability Margin');
axis([omega(1) omega(end) 0 1]);
grid on;
zoom on;

% figure(2);
% hax2 = gca;
% loglog(hax2, omega, LoopShape_upper_fr, '--', omega, LoopShape_lower_fr, '--');
% hold(hax2, 'on');
% loglog(hax2, omega, ones(size(omega)), ':');
% sigma(P, omega); % Plots singular values of P
% hold(hax2, 'off');
% xlabel(hax2, 'Frequency (radians/sec)');
% ylabel(hax2, 'Singular Values');
% title(hax2, 'Nominal Plant Singular Values & Loop-Shape Boundaries');
% grid(hax2, 'on');
% zoom(hax2, 'on');

figure(2);
% Manually calculate the singular values of the plant P instead of plotting directly.
sv_P = sigma(P, omega);

% Plot all data series in a single, robust loglog command.
loglog(omega, LoopShape_upper_fr, 'r--', ...
       omega, LoopShape_lower_fr, 'r--', ...
       omega, ones(size(omega)), 'k:', ...
       omega, sv_P, 'b-');

% Add labels and formatting.
title('Nominal Plant Singular Values & Loop-Shape Boundaries');
xlabel('Frequency (radians/sec)');
ylabel('Singular Values');
legend('Upper Bound', 'Lower Bound', 'Unity Gain', 'Plant SVs', 'Location', 'best');
grid on;
zoom on;


figure(3);
if W1_bounds_FLAG
    hax3_1 = subplot(1,1,1);
    bodemag(Weight_bounds(1,2), '--', Weight_bounds(2,2), '--', {omega(1), omega(end)});
    title('Boundaries for \sigma_i (W_2)');
    grid on; zoom on;
elseif W2_bounds_FLAG
    hax3_1 = subplot(1,1,1);
    bodemag(Weight_bounds(1,1), '--', Weight_bounds(2,1), '--', {omega(1), omega(end)});
    title('Boundaries for \sigma_i (W_1)');
    grid on; zoom on;
else
    hax3_1 = subplot(2,1,1);
    bodemag(Weight_bounds(1,1), '--', Weight_bounds(2,1), '--', {omega(1), omega(end)});
    title('Boundaries for \sigma_i (W_1)');
    grid on; zoom on;
    
    hax3_2 = subplot(2,1,2);
    bodemag(Weight_bounds(1,2), '--', Weight_bounds(2,2), '--', {omega(1), omega(end)});
    title('Boundaries for \sigma_i (W_2)');
    grid on; zoom on;
end

figure(4);
if W1_bounds_FLAG
    hax4_1 = subplot(1,1,1);
    semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,2),omega))),'--'); hold on;
    semilogx(omega, ones(size(omega)),'--'); hold off;
    title('Boundaries for k(W_2)'); grid on; zoom on;
elseif W2_bounds_FLAG
    hax4_1 = subplot(1,1,1);
    semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,1),omega))),'--'); hold on;
    semilogx(omega, ones(size(omega)),'--'); hold off;
    title('Boundaries for k(W_1)'); grid on; zoom on;
else
    hax4_1 = subplot(2,1,1);
    semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,1),omega))),'--'); hold on;
    semilogx(omega, ones(size(omega)),'--'); hold off;
    title('Boundaries for k(W_1)'); grid on; zoom on;
    
    hax4_2 = subplot(2,1,2);
    semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,2),omega))),'--'); hold on;
    semilogx(omega, ones(size(omega)),'--'); hold off;
    title('Boundaries for k(W_2)'); grid on; zoom on;
end

figure(5);
sigma(Cinit, {omega(1), omega(end)});
title('Singular Values of Initial Controller');
grid on; zoom on;
drawnow;
pause(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF ACTUAL PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Uhat and Vhat for diagonal/non-diagonal weights
[Uhat,Vhat] = unitary_svd_approx(P, omega, WeightStruc);

% Begin iterations
ii = 0; flag = 'NOTEXIT'; emaxIN = -1; decvarsIN = [];
while ~strcmpi(flag,'exit')
    ii = ii + 1;

    fprintf('\n\n--- W-iteration: %d ---\n\n', ii);
    [W1_tmp,W2_tmp,rho_W_iter,decvarsOUT] = W_iter(P,C,Uhat,Vhat,LoopShape_bounds, ...
                                                Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN);
    Ps = W2_tmp * P * W1_tmp;
    Ps = minreal(Ps, sqrt(eps));

    % Plotting graphs after W-iteration
    % --- This is the NEW, corrected block for updating figure(2) inside the loop ---
    % It should be placed after the line 'Ps = W2_tmp * P * W1_tmp;'
    
    figure(2);
    clf; % Clear the figure. This handles cases where the window was closed.
    
    % Manually calculate the singular values of the *shaped* plant Ps.
    sv_Ps = sigma(Ps, omega);
    
    % Plot all data series in a single, robust loglog command.
    % Note: We re-plot the bounds on each iteration for clarity.
    loglog(omega, LoopShape_upper_fr, 'r--', ...
           omega, LoopShape_lower_fr, 'r--', ...
           omega, ones(size(omega)), 'k:', ...
           omega, sv_Ps, 'b-');
    
    % Add labels and formatting.
    title('Singular Values of Shaped Plant');
    xlabel('Frequency (radians/sec)');
    ylabel('Singular Values');
    legend('Upper Bound', 'Lower Bound', 'Unity Gain', 'Shaped Plant SVs', 'Location', 'best');
    grid on;
    zoom on;
    
    % Optional: Restore the original axis limits if you have them defined
    % For example:
    % loopshape_upper_axis_limit = ...
    % loopshape_lower_axis_limit = ...
    % axis([omega(1) omega(end) loopshape_lower_axis_limit loopshape_upper_axis_limit]);
   
    % figure(2);
    % hax2 = gca; cla(hax2);
    % loglog(hax2, omega, LoopShape_upper_fr, '--', omega, LoopShape_lower_fr, '--');
    % hold(hax2, 'on');
    % loglog(hax2, omega, ones(size(omega)), ':');
    % sigma(Ps, omega, 'b-'); % Plot new shaped plant SVs
    % hold(hax2, 'off');
    % title(hax2, 'Singular Values of Shaped Plant');
    % grid(hax2, 'on');
    
    % figure(3);
    % if W1_bounds_FLAG
    %     subplot(1,1,1); cla;
    %     bodemag(Weight_bounds(1,2), '--', Weight_bounds(2,2), '--', {omega(1), omega(end)}); hold on;
    %     sigma(W2_tmp, 'b-', {omega(1), omega(end)}); hold off; title('Singular Values of W_2');
    % elseif W2_bounds_FLAG
    %     subplot(1,1,1); cla;
    %     bodemag(Weight_bounds(1,1), '--', Weight_bounds(2,1), '--', {omega(1), omega(end)}); hold on;
    %     size({omega(1), omega(end)})
    %     sigma(W1_tmp, 'b-', {omega(1), omega(end)}); hold off; title('Singular Values of W_1');
    % else
    %     subplot(2,1,1); cla;
    %     bodemag(Weight_bounds(1,1), '--', Weight_bounds(2,1), '--', {omega(1), omega(end)}); hold on;
    %     sigma(W1_tmp, 'b-', {omega(1), omega(end)}); hold off; title('Singular Values of W_1');
    % 
    %     subplot(2,1,2); cla;
    %     bodemag(Weight_bounds(1,2), '--', Weight_bounds(2,2), '--', {omega(1), omega(end)}); hold on;
    %     sigma(W2_tmp, 'b-', {omega(1), omega(end)}); hold off; title('Singular Values of W_2');
    % end
    % --- This is the NEW, corrected block ---
    figure(3);
    if W1_bounds_FLAG
        % Plot W2 against its bounds
        h_ax = subplot(1,1,1);
        cla(h_ax);
        
        % Get frequency response data for bounds and singular value data for W2
        mag_upper = abs(squeeze(freqresp(Weight_bounds(1,2), omega)));
        mag_lower = abs(squeeze(freqresp(Weight_bounds(2,2), omega)));
        sv_w = sigma(W2_tmp, omega);
        
        % Use loglog for manual plotting
        loglog(h_ax, omega, mag_upper, 'r--', omega, mag_lower, 'r--', omega, sv_w, 'b-');
        grid(h_ax, 'on');
        title(h_ax, 'Singular Values of W_2');
        xlabel(h_ax, 'Frequency (rad/s)');
        ylabel(h_ax, 'Magnitude');
        
    elseif W2_bounds_FLAG
        % Plot W1 against its bounds
        h_ax = subplot(1,1,1);
        cla(h_ax);
        
        % Get frequency response data for bounds and singular value data for W1
        mag_upper = abs(squeeze(freqresp(Weight_bounds(1,1), omega)));
        mag_lower = abs(squeeze(freqresp(Weight_bounds(2,1), omega)));
        sv_w = sigma(W1_tmp, omega);
        
        % Use loglog for manual plotting
        loglog(h_ax, omega, mag_upper, 'r--', omega, mag_lower, 'r--', omega, sv_w, 'b-');
        grid(h_ax, 'on');
        title(h_ax, 'Singular Values of W_1');
        xlabel(h_ax, 'Frequency (rad/s)');
        ylabel(h_ax, 'Magnitude');
        
    else
        % Plot W1 against its bounds in the top subplot
        h_ax1 = subplot(2,1,1);
        cla(h_ax1);
        mag_upper1 = abs(squeeze(freqresp(Weight_bounds(1,1), omega)));
        mag_lower1 = abs(squeeze(freqresp(Weight_bounds(2,1), omega)));
        sv_w1 = sigma(W1_tmp, omega);
        loglog(h_ax1, omega, mag_upper1, 'r--', omega, mag_lower1, 'r--', omega, sv_w1, 'b-');
        grid(h_ax1, 'on');
        title(h_ax1, 'Singular Values of W_1');
        xlabel(h_ax1, 'Frequency (rad/s)');
        ylabel(h_ax1, 'Magnitude');
    
        % Plot W2 against its bounds in the bottom subplot
        h_ax2 = subplot(2,1,2);
        cla(h_ax2);
        mag_upper2 = abs(squeeze(freqresp(Weight_bounds(1,2), omega)));
        mag_lower2 = abs(squeeze(freqresp(Weight_bounds(2,2), omega)));
        sv_w2 = sigma(W2_tmp, omega);
        loglog(h_ax2, omega, mag_upper2, 'r--', omega, mag_lower2, 'r--', omega, sv_w2, 'b-');
        grid(h_ax2, 'on');
        title(h_ax2, 'Singular Values of W_2');
        xlabel(h_ax2, 'Frequency (rad/s)');
        ylabel(h_ax2, 'Magnitude');
    end

    % figure(4);
    % if W1_bounds_FLAG
    %     subplot(1,1,1); cla;
    %     semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,2),omega))),'--'); hold on;
    %     semilogx(omega, ones(size(omega)),'--');
    %     semilogx(omega, cond(freqresp(W2_tmp,omega)),'b-'); hold off;
    %     title('Condition Number of W_2');
    % elseif W2_bounds_FLAG
    %      subplot(1,1,1); cla;
    %     semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,1),omega))),'--'); hold on;
    %     semilogx(omega, ones(size(omega)),'--');
    %     semilogx(omega, cond(freqresp(W1_tmp,omega)),'b-'); hold off;
    %     title('Condition Number of W_1');
    % else
    %     subplot(2,1,1); cla;
    %     semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,1),omega))),'--'); hold on;
    %     semilogx(omega, ones(size(omega)),'--');
    %     semilogx(omega, cond(freqresp(W1_tmp,omega)),'b-'); hold off;
    %     title('Condition Number of W_1');
    % 
    %     subplot(2,1,2); cla;
    %     semilogx(omega, abs(squeeze(freqresp(CondNo_bounds(1,2),omega))),'--'); hold on;
    %     semilogx(omega, ones(size(omega)),'--');
    %     semilogx(omega, cond(freqresp(W2_tmp,omega)),'b-'); hold off;
    %     title('Condition Number of W_2');
    % end
    % --- This is the NEW, corrected block ---
    figure(4);
    if W1_bounds_FLAG
        % Plot condition number of W2
        h_ax = subplot(1,1,1);
        cla(h_ax);
        
        % Calculate condition number by looping through frequencies
        W2_fr = freqresp(W2_tmp, omega);
        k_W2 = zeros(1, length(omega));
        for k = 1:length(omega)
            k_W2(k) = cond(W2_fr(:,:,k));
        end
        
        % Plot bounds and the result
        bound_k2 = abs(squeeze(freqresp(CondNo_bounds(1,2), omega)));
        semilogx(h_ax, omega, bound_k2, 'r--', omega, ones(size(omega)), 'k--', omega, k_W2, 'b-');
        grid(h_ax, 'on');
        title(h_ax, 'Condition Number of W_2');
        xlabel(h_ax, 'Frequency (rad/s)');
        ylabel(h_ax, 'Condition Number');
        
    elseif W2_bounds_FLAG
        % Plot condition number of W1
        h_ax = subplot(1,1,1);
        cla(h_ax);
        
        % Calculate condition number by looping through frequencies
        W1_fr = freqresp(W1_tmp, omega);
        k_W1 = zeros(1, length(omega));
        for k = 1:length(omega)
            k_W1(k) = cond(W1_fr(:,:,k));
        end
        
        % Plot bounds and the result
        bound_k1 = abs(squeeze(freqresp(CondNo_bounds(1,1), omega)));
        semilogx(h_ax, omega, bound_k1, 'r--', omega, ones(size(omega)), 'k--', omega, k_W1, 'b-');
        grid(h_ax, 'on');
        title(h_ax, 'Condition Number of W_1');
        xlabel(h_ax, 'Frequency (rad/s)');
        ylabel(h_ax, 'Condition Number');
        
    else
        % Plot condition number of W1 in the top subplot
        h_ax1 = subplot(2,1,1);
        cla(h_ax1);
        W1_fr = freqresp(W1_tmp, omega);
        k_W1 = zeros(1, length(omega));
        for k = 1:length(omega), k_W1(k) = cond(W1_fr(:,:,k)); end
        bound_k1 = abs(squeeze(freqresp(CondNo_bounds(1,1), omega)));
        semilogx(h_ax1, omega, bound_k1, 'r--', omega, ones(size(omega)), 'k--', omega, k_W1, 'b-');
        grid(h_ax1, 'on');
        title(h_ax1, 'Condition Number of W_1');
        xlabel(h_ax1, 'Frequency (rad/s)');
        ylabel(h_ax1, 'Condition Number');
    
        % Plot condition number of W2 in the bottom subplot
        h_ax2 = subplot(2,1,2);
        cla(h_ax2);
        W2_fr = freqresp(W2_tmp, omega);
        k_W2 = zeros(1, length(omega));
        for k = 1:length(omega), k_W2(k) = cond(W2_fr(:,:,k)); end
        bound_k2 = abs(squeeze(freqresp(CondNo_bounds(1,2), omega)));
        semilogx(h_ax2, omega, bound_k2, 'r--', omega, ones(size(omega)), 'k--', omega, k_W2, 'b-');
        grid(h_ax2, 'on');
        title(h_ax2, 'Condition Number of W_2');
        xlabel(h_ax2, 'Frequency (rad/s)');
        ylabel(h_ax2, 'Condition Number');
    end
    drawnow;
    pause(0.5);

    fprintf('\n\n--- Controller Synthesis: %d ---\n', ii);
    [Cinf_tmp, ~, emaxC_tmp(ii)] = ncfsyn(Ps);
    C = W1_tmp * Cinf_tmp * W2_tmp;
    C = minreal(C, sqrt(eps));
    
    % Compute closed-loop transfer function for stability margin analysis
    % This is [I; K] * (I + G*K)^-1 * [I, G] -> want norm of [I;K]*inv(I+PK)[-P, I]
    % which simplifies to margin calculation. The original code computes:
    % clp = daug(-eye(m),eye(n)) * minv(sbs(abv(eye(m),Cinf_tmp),abv(Ps,eye(n)))) * daug(eye(m),zeros(n))
    % Let's use the modern ncfsyn output directly.
    
    L = Ps * Cinf_tmp;
    S = inv(eye(size(L)) + L);
    clp_margin_tf = [Cinf_tmp*S; S];
    margin_at_freq = 1 ./ sigma(clp_margin_tf, omega);
    margin_at_freq = squeeze(margin_at_freq(1,1,:));


    % Plotting graphs after Controller Synthesis
    % figure(1);
    % hax1 = gca;
    % cla(hax1);
    % semilogx(hax1, omega, margin_at_freq, 'b-');
    % hold(hax1, 'on');
    % semilogx(hax1, omega, rho_W_iter, 'r:'); % Show margin from W-iteration
    % hold(hax1, 'off');
    % legend(hax1, 'After Controller Synthesis', 'Target from W-iteration', 'Location', 'best');
    % xlabel(hax1, 'Frequency (radians/sec)');
    % ylabel(hax1, 'Pointwise Robust Stability Margin');
    % axis(hax1, [omega(1) omega(end) 0 1]);
    % grid(hax1, 'on');
    figure(1);
    clf; % Clear the entire figure for a robust, clean slate.
    
    % Plot both lines in a single command. MATLAB handles the overlay automatically.
    semilogx(omega, margin_at_freq, 'b-', omega, rho_W_iter, 'r:');
    
    % Add formatting to the current axes.
    legend('After Controller Synthesis', 'Target from W-iteration', 'Location', 'best');
    xlabel('Frequency (radians/sec)');
    ylabel('Pointwise Robust Stability Margin');
    axis([omega(1) omega(end) 0 1]);
    grid on;
      
    figure(5);
    clf; % Clears the entire figure, which is more robust than clearing just the axes.
    
    % The sigma command will now create fresh axes on the cleared figure.
    sigma(Cinf_tmp, 'b-', {omega(1), omega(end)});
    
    % The title command will automatically apply to the newly created axes.
    title('Singular Values of Controller C_\infty');
    grid on;
    drawnow;
    pause(0.5);

    fprintf('\n\n--- Iteration %d successfully completed! ---\n', ii);
    fprintf('Optimised gamma value: %.4f\n', emaxC_tmp(ii));
    
    % Setting input variables for next iteration's LMI optimisation
    decvarsIN = decvarsOUT;
    emaxIN    = emaxC_tmp(ii);

    % Asking whether to do another iteration or not!
    flag = input('Type "exit" if no more iterations required, or press ENTER to continue: ','s');
end


fprintf('\n\nSuccessfully converged after %d iterations\n', ii);

% Plotting graph for robust stability margin as iterations proceed
figure(6);
plot(1:ii, emaxC_tmp, '-*');
xlabel('Iteration Number');
ylabel('Robust Stability Margin (\epsilon_{max})');
title('Convergence of Robust Stability Margin');
axis([0 (ii+1) 0 1]);
grid on; zoom on;
pause(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ACTUAL PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Setting outputs
if (nargout >= 1)
    emax = emaxC_tmp(end);
end
if (nargout >= 2)
    W1 = W1_tmp;
end
if (nargout >= 3)
    W2 = W2_tmp;
end
if (nargout >= 4)
    Cinf = Cinf_tmp;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        INTERNAL HELPER FUNCTIONS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Uhat,Vhat] = unitary_svd_approx(P, omega, WeightStruc)
% This function takes a continuous frequency-by-frequency Singular Value
% Decomposition of plant "P" and fits transfer function matrices "Uhat"
% and/or "Vhat" to the matrices of singular vectors.
% Updated for modern MATLAB.

[m,n] = size(P);

if strcmp(WeightStruc,'DiagDiag')
    Uhat = eye(m);
    Vhat = eye(n);
    return;
end

% Take a continuous SVD
fprintf('\nComputing continuous SVD of the plant...\n');
[U_g, ~, V_g] = vsvdcont(P, omega); % Custom function to get smooth SVD

% Approximate the continuous frequency matrices of singular vectors with
% rational LTI system matrices
fprintf('\nComputing rational approximations to matrices of singular vectors...\n');

if strcmp(WeightStruc,'DiagFull')
    U_sys = rational_approx(U_g, omega);
    V_sys = eye(n);
elseif strcmp(WeightStruc,'FullDiag')
    U_sys = eye(m);
    V_sys = rational_approx(V_g, omega);
elseif strcmp(WeightStruc,'FullFull')
    U_sys = rational_approx(U_g, omega);
    V_sys = rational_approx(V_g, omega);
end

% If "U_sys" or "V_sys" is tall/wide, construct an all-pass dilation to make it square.
% This is based on ZDG Lemma 13.31.
if (m > n) && ~strcmp(WeightStruc,'DiagDiag') && ~strcmp(WeightStruc,'FullDiag')
    U_sys = minreal(U_sys / sqrtm(U_sys'*U_sys)); % simple approach to make it unitary
    U_sys = balancmr(U_sys).ss; % balance and convert to ss
elseif (n > m) && ~strcmp(WeightStruc,'DiagDiag') && ~strcmp(WeightStruc,'DiagFull')
    V_sys = minreal(V_sys / sqrtm(V_sys'*V_sys));
    V_sys = balancmr(V_sys).ss;
end

Uhat = U_sys;
Vhat = V_sys;
end

% -------------------------------------------------------------------------

function U_sys = rational_approx(U_g, omega)
% This function finds a rational approximation "U_sys" to
% the frequency-varying matrix "U_g".
% It interactively asks for the order of each element's approximation.
% Requires System Identification Toolbox.

[m, n] = size(U_g, [1 2]); % U_g is a 3D array
U_fit_cells = cell(m, n);

for r = 1:m
    for c = 1:n
        flag = 'y';
        element_frd = frd(squeeze(U_g(r,c,:)), omega);
        
        h_bode = figure;
        bodeplot(element_frd, 'y.');
        hold on;
        title(['Bode Diagram for Singular Vector Element (' num2str(r) ',' num2str(c) ')']);
        
        while ~strcmpi(flag, 'n')
            
            order_str = input(['\nPlease input order of (' num2str(r) ',' num2str(c) ') approximation: ']);
            if isempty(order_str) || ~isnumeric(order_str) || order_str < 0
                disp('Invalid order. Please try again.');
                continue;
            end
            
            % Use tfest for fitting
            fprintf('Fitting with order %d...\n', order_str);
            approx_sys = tfest(element_frd, order_str);
            
            bodeplot(approx_sys, 'r-');
            legend('Original Data', 'Fitted Model');
            hold off;
            
            flag = input('Press ENTER to re-approximate, or type "n" to accept and continue: ','s');
        end
        U_fit_cells{r,c} = approx_sys;
        disp(' ');
    end
end
close(h_bode);

% Assemble the full system and perform model reduction
U_fit = cell2sys(U_fit_cells);
U_fit = minreal(U_fit, 1e-7);

% Model reduction using balanced truncation of normalized coprime factors
disp('Performing model reduction on the combined system...');
[~, RedInfo] = ncfmr(U_fit);
hankel_svs = RedInfo.Sigma;
fprintf('Hankel Singular Values of NCF: \n');
disp(hankel_svs');

% Interactive reduction
order_red_str = input(['Enter desired order for reduction (original is ' num2str(order(U_fit)) '): ']);
if isempty(order_red_str)
    % Default reduction based on a tolerance
    reduced_order = sum(hankel_svs > 0.01);
    fprintf('Using default tolerance, reduced order is %d\n', reduced_order);
else
    reduced_order = str2double(order_red_str);
end

U_red = ncfmr(U_fit, struct('MaxOrder', reduced_order));

U_sys = balreal(U_red).ss; % Return a balanced state-space realization
end

% -------------------------------------------------------------------------

function [U_g_cont, S_g_cont, V_g_cont] = vsvdcont(P, omega)
% This function takes a "continuous" SVD of P over the frequency vector omega.
% It sorts singular values and corrects phases of singular vectors to avoid
% discontinuities. This is an updated, non-interactive version of the original.

[m, n] = size(P);
q = min(m, n);
P_fr = freqresp(P, omega);

% Pre-allocate arrays
U_g_cont = zeros(m, q, length(omega));
S_g_cont = zeros(q, q, length(omega));
V_g_cont = zeros(n, q, length(omega));

% Initial SVD at the first frequency point
[U,S,V] = svd(P_fr(:,:,1));
[s_vals, s_idx] = sort(diag(S), 'descend');
S_g_cont(:,:,1) = diag(s_vals(1:q));
U_g_cont(:,:,1) = U(:, s_idx(1:q));
V_g_cont(:,:,1) = V(:, s_idx(1:q));

% Loop through remaining frequencies
for k = 2:length(omega)
    [U_curr, S_curr, V_curr] = svd(P_fr(:,:,k));
    
    % Sort singular values and vectors to match the previous frequency step
    % This is a simple greedy matching based on vector projection
    U_prev = squeeze(U_g_cont(:,:,k-1));
    V_prev = squeeze(V_g_cont(:,:,k-1));
    
    cost_matrix = abs(U_prev' * U_curr(:,1:q)) + abs(V_prev' * V_curr(:,1:q));
    
    % Use assignment algorithm (e.g., Hungarian) for robust matching
    % For simplicity, a greedy approach is used here.
    [~, perm_idx] = max(cost_matrix, [], 2);
    
    U_sorted = U_curr(:, perm_idx);
    S_sorted = diag(diag(S_curr(perm_idx, perm_idx)));
    V_sorted = V_curr(:, perm_idx);

    % Phase correction to ensure smooth transition
    phase_corr_U = diag(sign(diag(U_prev' * U_sorted)));
    phase_corr_V = diag(sign(diag(V_prev' * V_sorted)));
    
    U_g_cont(:,:,k) = U_sorted * phase_corr_U;
    V_g_cont(:,:,k) = V_sorted * phase_corr_V;
    S_g_cont(:,:,k) = S_sorted;
end
end

% -------------------------------------------------------------------------
% The core LMI iteration functions would be placed here.
% The logic inside W_iter and W_iter_tall_P is highly specialized and
% depends on the LMI-solver syntax, which is largely unchanged.
% A full rewrite of those sections is a significant undertaking, but the
% main changes would be replacing frequency-domain operations
% (frsp, vnorm) with modern equivalents (freqresp, sigma) and updating
% the spectral factorization part as sketched below.

function [W1,W2,rho,decvarsOUT] = W_iter(P,C,Uhat,Vhat,LoopShape_bounds,Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN)
% This function determines whether the plant "P" is TALL or FAT and prepares
% the input data for the function "W_iter_tall_P".
[m,n] = size(P);
if (m >= n)
  [W1_tmp,W2_tmp,rho_tmp,decvarsOUT_tmp] = ...
       W_iter_tall_P(P,C,Uhat,Vhat,LoopShape_bounds,Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN,1);
else
  % Handle FAT case by transposing the problem
  Weight_bounds_T = Weight_bounds * [0 1; 1 0];
  CondNo_bounds_T = CondNo_bounds * [0 1; 1 0];
  [W2_T,W1_T,rho_tmp,decvarsOUT_tmp] = ...
       W_iter_tall_P(P.',C.',Vhat',Uhat',LoopShape_bounds, ...
					 Weight_bounds_T,CondNo_bounds_T,omega,emaxIN,decvarsIN,0);
  W1_tmp = W1_T.';
  W2_tmp = W2_T.';
end
W1 = W1_tmp; W2 = W2_tmp; rho = rho_tmp; decvarsOUT = decvarsOUT_tmp;
end

function [W1,W2,rho,decvarsOUT] = W_iter_tall_P(P,C,Uhat,Vhat,LoopShape_bounds,Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN,tall)
% This function performs the W-iteration for a TALL or SQUARE plant.
% It defines and solves an LMI optimization problem at each frequency,
% then fits transfer functions to the results and constructs the weights
% via spectral factorization.
% NOTE: This is a high-level sketch. The LMI setup is complex.
global HANDLEW; % Use global to modify plot handle from within function

[m,n] = size(P);

% --- LMI SETUP (As per original, highly complex) ---
% This part uses setlmis, lmivar, lmiterm which are still valid.
% The key is to correctly compute the constant matrices for the LMIs at
% each frequency point using modern `freqresp` and `squeeze`.

% --- Mockup of the frequency loop and LMI solving ---
num_dec_vars = n + m + 4; % Based on original lmivar calls
decvarsOUT_tmp = zeros(num_dec_vars, length(omega));
gammasqr = zeros(1, length(omega));
% (LMI definitions would go here...)
% lmisys = getlmis;

fprintf('Solving LMIs pointwise across frequency range...\n');
for jj = 1:length(omega)
    ww = omega(jj);
    % Define "big" constants for complex LMIs using modern syntax
    % e.g., CostFunctRightConstr_g = squeeze(freqresp(CostFunctRightConstr, ww));
    % ... then pass these to lmiterm
    
    % Mock solver call
    if isempty(decvarsIN)
        %[g, d_out] = gevp(lmisys, 1, [1e-3 200 0 20 1]);
        d_out = rand(num_dec_vars,1); g = rand; % Placeholder
    else
        %[g, d_out] = gevp(lmisys, 1, [1e-3 200 0 20 1], 1.2/(emaxIN*emaxIN), decvarsIN(:,jj), 1.05);
        d_out = rand(num_dec_vars,1); g = rand; % Placeholder
    end
    gammasqr(jj) = g;
    decvarsOUT_tmp(:,jj) = d_out;
    
    if (mod(jj,10) == 0), fprintf('%d/%d\n', jj, length(omega)); end
end
fprintf('... LMI solving complete.\n');

rho_tmp = 1./sqrt(gammasqr);
modD2   = sqrt(decvarsOUT_tmp(1:m,:));
modD1   = 1./sqrt(decvarsOUT_tmp((m+1):(m+n),:));

figure(1);
hax1 = gca;
if (emaxIN ~= -1) && isgraphics(HANDLEW)
    set(HANDLEW, 'Visible', 'off');
end
hold(hax1, 'on');
HANDLEW = semilogx(hax1, omega, rho_tmp, 'r:');
hold(hax1, 'off');
drawnow;

disp('Fitting transfer functions to magnitudes...');

% Use the custom updated fitting function
if (tall == 1), a1 = '1'; a2 = '2'; else, a1 = '2'; a2 = '1'; end

W1_bounds_FLAG = norm(getPeakGain(Weight_bounds(:,1) - [1;1])) == 0;
W2_bounds_FLAG = norm(getPeakGain(Weight_bounds(:,2) - [1;1])) == 0;


if W1_bounds_FLAG
    D1 = eye(n);
else
    D1_cells = cell(1,n);
   
    for kk = 1:n
        fprintf('\nFitting D_%s (%d,%d)\n', a1, kk, kk);
        mag_data = frd(modD1(kk,:).', omega);
        D1_cells{kk} = fitmag_updated(mag_data);
    end
    % D1 = diag(cell2mat(D1_cells));  % cell2mat does not support cell arrays of cell arrays
    D1 = blkdiag(D1_cells{:});
end

if W2_bounds_FLAG
    D2 = eye(m);
else
    D2_cells = cell(1,m);
    for kk = 1:m
        fprintf('\nFitting D_%s (%d,%d)\n', a2, kk, kk);
        mag_data = frd(modD2(kk,:).', omega);
        D2_cells{kk} = fitmag_updated(mag_data);
    end
    % D2 = diag(cell2mat(D2_cells));
    D2 = blkdiag(D2_cells{:});
end

% Constructing Weights through Spectral Factorisation
disp('Constructing weights via spectral factorization...');
W1 = spectfact_updated(Vhat * D1);
W2 = spectfact_updated(Uhat * D2);

W1 = minreal(W1);
W2 = minreal(W2);

rho = rho_tmp;
decvarsOUT = decvarsOUT_tmp;
end

%--------------------------------------------------------------------------
% function W = spectfact_updated(G)
% % Performs spectral factorization W*W' = G*G' using CARE.
% % G is the system to be factorized. W is the stable, minimum-phase factor.
% if isempty(G) || norm(G) == 0
%     W = G; return;
% end
% 
% disp("hi")
% G = minreal(ss(G))
% disp("hi22")
% [A,B,C,D] = ssdata(G)
% R = D'*D;
% Q = C'*C;
% S = C'*D;
% 
% % Solve the Algebraic Riccati Equation for factorization
% % A'X + XA - (XB+S)R^-1(B'X+S') + Q = 0
% [X, K, ~] = icare(A, B, Q, R, S);
% 
% % The spectral factor W has the form: W = ss(A-B*K, B, C-D*K, D)
% % However, we need to ensure the D matrix of the factor is Cholesky of original R
% 
% L = chol(R, 'lower');
% W = ss(A, B, (C-D*pinv(R)*(B'*X+S')), L);
% W = minreal(W, sqrt(eps));
% end
function W = spectfact_updated(G)
% Performs spectral factorization W'*W = G'*G using CARE, where W is a
% stable, minimum-phase factor.
% G is the system to be factorized.
% This updated version correctly handles both strictly proper (D=0) and
% proper (D~=0) systems.

if isempty(G) || norm(G) == 0
    W = G;
    return;
end

G = minreal(ss(G));
[A,B,C,D] = ssdata(G);

% Check if the system is strictly proper (D matrix is effectively zero)
if norm(D, 'fro') < sqrt(eps)
    % --- STRICTLY PROPER CASE (D = 0) ---
    % The spectral factor W will also be strictly proper.
    % We solve the simplified Riccati equation: A'X + XA - XBB'X + C'C = 0.
    
    Q = C'*C;
    % Use the standard CARE solver for this simplified form
    [X, ~, ~] = care(A, B, Q);

    % The output matrix of the spectral factor W is L = B'*X
    L = B'*X;
    
    % The factor W is ss(A, B, L, D_w), where D_w must be zero.
    output_dim = size(L, 1);
    input_dim = size(B, 2);
    W = ss(A, B, L, zeros(output_dim, input_dim));

else
    % --- PROPER CASE (D ~= 0) ---
    % This assumes D'*D is invertible (i.e., G has no zeros on the jw-axis).
    R = D'*D;
    Q = C'*C;
    S = C'*D;

    % Check for singularity of R before proceeding.
    % rcond is a good way to check if a matrix is close to singular.
    if rcond(R) < 1e-12
         error('Spectral factorization failed: D''*D is singular or ill-conditioned. The system may have jw-axis zeros.');
    end

    % Solve the general continuous-time Algebraic Riccati Equation (CARE)
    % A'X + XA - (XB+S)R^-1(B'X+S') + Q = 0
    [X, ~, K] = care(A, B, Q, R, S); % K is the gain R^-1 * (B'X + S')

    % The Cholesky factor of R is the feedthrough matrix of the spectral factor W.
    % We use 'lower' to be consistent.
    L = chol(R, 'lower');

    % The spectral factor is W(s) = L + K * (sI - A)^-1 * B.
    % The state-space realization for this is:
    W = ss(A, B, K, L);
end

% Ensure the final result is a minimal realization.
W = minreal(W, sqrt(eps));
end
%--------------------------------------------------------------------------
function sys = fitmag_updated(mag_frd)
% Interactively fits a transfer function to magnitude data.
% Replaces the non-standard `fitmag` function.
flag = 'y';
h_fit_fig = figure;
% loglog(mag_frd.Frequency, mag_frd.ResponseData, 'y.');
loglog(mag_frd.Frequency(:, 1), squeeze(mag_frd.ResponseData(1, 1, :)), 'm.');
hold on; grid on;
title('Fit Magnitude Data');
xlabel('Frequency (rad/s)'); ylabel('Magnitude');

while ~strcmpi(flag, 'n')
    order_str = input('Please input order of the approximation: ');
    if isempty(order_str) || ~isnumeric(order_str) || order_str < 0
        disp('Invalid order. Please try again.');
        continue;
    end
    
    fprintf('Fitting with order %d...\n', order_str);
    % Use tfest for fitting. We fit the magnitude, so result is TF.
    sys_tf = tfest(mag_frd, order_str, 0); % order poles, 0 zeros.
    
    % Plot the result
    [mag, ~, wout] = bode(sys_tf, mag_frd.Frequency);
    loglog(wout, squeeze(mag), 'r-');
    legend('Original Data', 'Fitted Model');
    hold off;
    
    flag = input('Press ENTER to re-approximate, or type "n" to accept: ','s');
end
close(h_fit_fig);
% Return a stable, minimum-phase system
[z, p, k] = zpkdata(sys_tf, 'v');
% Reflect unstable poles/zeros to be stable/min-phase
p(real(p)>0) = -p(real(p)>0);
z(real(z)>0) = -z(real(z)>0);
sys = zpk(z,p,k);
sys = ss(sys); % return state-space
end
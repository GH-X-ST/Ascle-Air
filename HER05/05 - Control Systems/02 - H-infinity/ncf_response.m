function [K, CL, gamma] = ncf_response(P, W1, W2)
% Ingests pre-compensator and post-compensator weighting functions,
% produces a H-infinity controller via ncfsyn and returns the gamma value
% along with singular value and frequency response Bode plots.
    arguments   
        % required inputs
        P (1, 1) lti

        % optional arguments
        W1 (1, 1) lti = ss(1)
        W2 (1, 1) lti = ss(1)
        
    end
    
    % Display robustness characteristics before
    pre_margins = allmargin(W1 * P * W2);
    disp("Shaped plant margin: ")
    disp(pre_margins)

    % Synthesise controller using ncfsyn
    [K, CL, gamma] = ncfsyn(-P, W1, W2);  % ncfsyn assumes +ve feedback
    ncf_margins = allmargin(P*K);
    disp("NCFSYN margins: ")
    disp(ncf_margins)

    disp(['NCFSYN gamma value: ', num2str(gamma)]);

    % Plot sensitivity function and save as png
    figure;
    sigma(P, 'g', P * K, 'b--', W1 * P * W2, 'r:', {1e-3, 1e2});
    title("Sensitivity function")
    legend("Open loop", "NCFSYN design", "Compensator design");
    print("sigma_sens.png", "-dpng");

    % Plot system response
    figure;
    bodemag(P, 'g', P * K, 'b--', W1 * P * W2, 'r:', {1e-3, 1e2});
    legend("Open loop", "NCFSYN design", "Compensator design");
    title("System response")
    print("sys_resp.png", "-dpng");
    
     % Plot impulse response
    figure;
    impulse(feedback(P, K), 'b', P, 'r:', W1 * P * W2, 'g--', 15);
    title("Impulse response")
    ylim([-2, 2]);
    legend("NCFSYN", "Open loop", "Compensator design")
    print("impulse.png", "-dpng");

    % Plot step response
    figure;
    step(feedback(P, K), 'b', P, 'r:', W1 * P * W2, 'g--', 15);
    ylim([-2, 2]);
    title("Step response")
    legend("NCFSYN", "Open loop", "Compensator design")
    print("step_resp.png", "-dpng");

end
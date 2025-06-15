function [] = get_response_graphs(P, K)
% Ingests pre-compensator and post-compensator weighting functions,
% produces a H-infinity controller via ncfsyn and returns the gamma value
% along with singular value and frequency response Bode plots.
   
    
    % Display robustness characteristics 
    pre_margins = allmargin(P * K);
    disp(pre_margins)

    % Synthesise controller using ncfsyn


    % Plot sensitivity function and save as png
    figure;
    sigma(P, 'g', P * K, 'b--', {1e-4 1e3});
    title("Sensitivity function")
    yline(0, '--')
    legend("Open loop", "H-infinity design", '0 dB');
    print("sigma_sens_mix.png", "-dpng");

    % Plot system response
    figure;
    bode(P, 'g', P * K, 'b--', {1e-4, 1e3});
    yline(-180, '--')
    legend("Open loop", "H-infinity design", "-180 degrees");
    print("sys_resp_mix.png", "-dpng");
    
     % Plot impulse response
    figure;
    impulse(feedback(P, K), 'b', P, 'r:', 15);
    title("Impulse response")
    ylim([-2, 2]);
    legend("H-infinity design", "Open loop")
    print("impulse_mix.png", "-dpng");

    % Plot step response
    figure;
    step(feedback(P, K), 'b', P, 'r:', 15);
    ylim([-2, 2]);
    title("Step response")
    legend("H-infinity design", "Open loop")
    print("step_resp_mix.png", "-dpng");
    
    % Get Nyquist loop
    figure;
    np = nyquistplot(P*K);


end
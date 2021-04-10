function val=inputC(t)
    % Gaussian Pulse
    sg = 0.03;
    mu = 0.06^2;
    val = exp(-(t-mu)^2 / (2*sg^2));
end 

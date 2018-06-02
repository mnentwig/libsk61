function sk61example_td_estFundPeriod()
    % === load library into variable ===
    addpath('..'); % location to libsk61_0v1.m

    % optionally: turn it into global variable, then use the same global statement in functions
    global libsk61;
    libsk61 = libsk61_0v1();
    
    % === example function spanning [0..1[ range ===
    x = 0:10000;
    T = 123.456;    
    y = cos(x/T*2*pi);
    tEst_samples = libsk61.td.estFundPeriod(y, 1.99*T);
    fprintf(stdout, 'True period: %1.5d estimated: %1.5f\n', T, tEst_samples);
end
function sk61example_spline()
    close all;
    graphics_toolkit('gnuplot');
    
    % === load library into variable ===
    addpath('..'); % location to libsk61_0v1.m

    % optionally: turn it into global variable, then use the same global statement in functions
    global libsk61;
    libsk61 = libsk61_0v1();
    
    % === example function spanning [0..1[ range ===
    nSeg = 16;
    x = linspace(0, 1, 5001); x = x(1:end-1);
    y = cos(x*2*pi);
    
    % === sweep the number of segments and show the error ===
    fprintf(stdout, 'nSeg\tapproximation error (dB)\n');
    for nSeg = 2:16
        % the spline functions assume an input range [0..nSeg[
        xx = x * nSeg;
                
        c = libsk61.spline.create(xx, y, nSeg, [], []);
        ye = libsk61.spline.eval(xx, c);

        err_dB = 20*log10(norm(ye-y) / norm(y));
        fprintf(stdout, '%i\t%1.2f\n', nSeg, err_dB);
    end

    % === plot the last spline ===
    figure(); 
    subplot(2, 1, 1); hold on; leg = {};
    plot(xx, y, 'k', 'lineWidth', 3); leg{end+1} = 'original function';
    plot(xx, ye, 'r'); leg{end+1} = 'spline fit';
    legend(leg);    
    subplot(2, 1, 2); leg = {};
    plot(xx, 20*log10(abs(ye-y) + 1e-3)); leg{end+1} = 'error';
    ylabel('dB');
    legend(leg);
    
    % === example function with discontinuities ===
    y = y .* linspace(1, 1.5, numel(x));    
    mask = (xx >= 4) & (xx < 6);
    y(mask) = y(mask) * 2;
    
    figure(); hold on; leg = {};
    plot(xx, y, 'k', 'lineWidth', 3); leg{end+1} = 'original function';

    c = libsk61.spline.create(xx, y, nSeg, [], []);
    plot(xx, libsk61.spline.eval(xx, c)); leg{end+1} = 'cyclic spline';

    discontVal = [0, 6];
    discontSlope = [0, 4, 6];
    c = libsk61.spline.create(xx, y, nSeg, discontVal, discontSlope);
    plot(xx, libsk61.spline.eval(xx, c)); leg{end+1} = 'discontinuities set';
    
    legend(leg);    
end

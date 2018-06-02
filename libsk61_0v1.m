function lib = libsk61_0v1()
% usage: 
% lib = libsk61();
% use functions with lib.functionname, e.g. lib.spline.create(...)
% omit arguments to show usage information, e.g. lib.spline.create()
%
% lib contains functions:
%     spline.create
%     spline.eval
%     td2td.noncausalResampler
%     td2td.delay
%     td2td.dcBlock
%     td.estFundPeriod
    persistent self;
    if ~isempty(self) 
        lib = self;
        return; 
    end
    self = struct();
    self.spline = struct(...
        'create', @spline_create, ...
        'eval', @spline_eval);
    self.td2td = struct(...
        'noncausalResampler', @td2td_noncausalResampler, ...
        'delay', @td2td_delay, ...
        'dcBlock', @td2td_dcBlock);
    self.td = struct(...
        'estFundPeriod', @td_estFundPeriod);
    lib = self;
end

function y = spline_eval(x, c)
    if nargin == 0
        d();
        d('y = spline_eval(x, c)');
        d();
        d('    evaluates a 3rd order spline');
        d('x:');
        d('    location of points');
        d('    0 <= x < nSeg with nSeg = numel(c)/4');
        d('c:');
        d('    spline coefficients, see spline_create');
        return;
    end
    if nargin ~= 2
        error('need two arguments x and c');
    end
    
    intPart = floor(x); intPart = intPart(:);
    fracPart = (mod(x, 1) - 0.5)*2; fracPart = fracPart(:);
    if intPart < 0 error('out of range: x is below spline range'); end
    if 4*intPart+4 > numel(c) error('out of range: x is above spline range'); end
    if false
        figure(); hold on;
        h = plot(intPart, 'k'); set(h, 'lineWidth', 3);
        h = plot(fracPart, 'b'); set(h, 'lineWidth', 3)
        legend('integer part (bank index', 'fractional part (indep. variable)');;    
    end
    
    ixc3 = 4*intPart+1;
    ixc2 = 4*intPart+2;
    ixc1 = 4*intPart+3;
    ixc0 = 4*intPart+4;
    
    % evaluate via Horner scheme
    % z = d(ixc0) + d(ixc1) .* fracPart + d(ixc2) .* fracPart.^2 + d(ixc3) .* fracPart.^3;
    y = c(ixc3);
    y = y .* fracPart;
    y = y + c(ixc2);
    y = y .* fracPart;
    y = y + c(ixc1);
    y = y .* fracPart;
    y = y + c(ixc0);

    if size(x) ~= size(y)
        y = y .';
        assert(size(x) == size(y));
    end
end

function c = spline_create(x, y, nSeg, xDiscontLvl, xDiscontSlope)
    if nargin == 0
        d('');
        d('c = spline_create(x, y, n, xDiscontLvl, xDiscontSlope)');
        d('    construct a 3rd order spline approximating data');
        d('');
        d('nSeg:');
        d('    number of segments (e.g. nSeg=4 results in 4 segments, 5 support points)');
        d('c:');
        d('    coefficients of polynomial segments [c10 c11 c12 c13 c20 c21 c22 c23 c30');
        d('    where cik is xLocal^k in segment i and -1 <= xLocal < 1');
        d('x:');
        d('    support points of data');
        d('    0 <= x <= nSeg');
        d('y:');
        d('    data values');
        d('xDiscontLvl: default [1]');
        d('    x values (integer) where the function value may be discontinuous');
        d('    an empty list fits a cyclic spline');
        d('xDiscontSlope: default xDiscontLvl');
        d('    x values (integer) where the slope may be discontinuous');
        d('    an empty list fits a cyclic spline');
        return;
    end
    if nargin < 3 error('need at least x, y, n'); end
    if nargin < 4 xDiscontLvl = [1]; end
    if nargin < 5 xDiscontSlope = xDiscontLvl; end
    if floor(xDiscontLvl) ~= xDiscontLvl error('xDiscontLvl must be integer (can split only between poly segments)'); end
    if floor(xDiscontSlope) ~= xDiscontSlope error('xDiscontSlope must be integer (can split only between poly segments)'); end
    xDiscontLvl = mod(xDiscontLvl, nSeg);
    xDiscontSlope = mod(xDiscontSlope, nSeg);
    
    % *****************************************************
    % basis polynomials for 3rd order spline
    % *****************************************************
    p1=[1/4,0,-3/4,1/2];    % [p(-1), p(1), dp/dx(-1), dp/dx(1)] = [1 0 0 0]
    p2=[-1/4,0,3/4,1/2];    % [p(-1), p(1), dp/dx(-1), dp/dx(1)] = [0 1 0 0]
    p3=[1/4,-1/4,-1/4,1/4]; % [p(-1), p(1), dp/dx(-1), dp/dx(1)] = [0 0 1 0]
    p4=[1/4,1/4,-1/4,-1/4]; % [p(-1), p(1), dp/dx(-1), dp/dx(1)] = [0 0 0 1]    
    
    M = zeros(numel(x), 1);
    MM = [];
    nSeg = nSeg;
    ixc = 0;
    for ixSeg = 1 : nSeg
        splitLevelFlag = max(xDiscontLvl == ixSeg-1);
        splitSlopeFlag = max(xDiscontSlope == ixSeg-1);

        % === forwards looking segment ===
        ixp = 4*ixSeg-3:4*ixSeg;
        xEv = interp1([ixSeg-1:ixSeg], [-1, 1], x+1e-15, 'linear', 'extrap');
        mask = (xEv >= -1) & (xEv < 1);
        
        ixc = ixc + 1;
        ixCoeffLevelForw = ixc;
        M(mask, ixc) = polyval(p1, xEv(mask));
        MM(ixc, ixp) = p1;
        
        ixc = ixc + 1;
        ixCoeffSlopeForw = ixc;
        M(mask, ixc) = polyval(p3, xEv(mask));
        MM(ixc, ixp) = p3;
        
        % === backwards looking segment ===
        ixSegBack = mod(ixSeg-2, nSeg) + 1;
        ixp = 4*ixSegBack-3:4*ixSegBack;
        xEv = interp1([ixSegBack-1:ixSegBack], [-1, 1], x+1e-15, 'linear', 'extrap');    
        mask = (xEv >= -1) & (xEv < 1);
        
        if splitLevelFlag
            ixc = ixc + 1;
            ixcBack = ixc;
        else
            ixcBack = ixCoeffLevelForw;
        end
        M(mask, ixcBack) = polyval(p2, xEv(mask));
        MM(ixcBack, ixp) = p2;
        
        if splitSlopeFlag
            ixc = ixc + 1;
            ixcBack = ixc;
        else
            ixcBack = ixCoeffSlopeForw;
        end
        
        M(mask, ixcBack) = polyval(p4, xEv(mask));
        MM(ixcBack, ixp) = p4;
    end
    if size(M) ~= [numel(x), ixc]
        % note: apparently unnecessary - M gets resized above, even if mask has no true entry
        M(numel(x), ixc) = 0;
    end

    % *****************************************************
    % === calculate basis function coefficients ===
    % *****************************************************
    c = pinv(M) * y(:);
    
    % *****************************************************
    % === convert basis function coefficients to polynomial coefficients ===
    % *****************************************************
    c = MM.'*c;
end

function d(s)
    if nargin < 1 s = ''; end
    disp(['    ', s]);
end

function td = td2td_noncausalResampler(td, n)

    if nargin == 0
        d();
        d('td = td2td_noncausalResampler(td, n)');
        d();
        d('    resamples a cyclic signal by ideal lowpass filtering');
        d('td:');
        d('    input signal');
        d('n: (integer)');
        d('    number of samples in output');
        return;
    end

    assert(n > 0);
    nOld = numel(td); assert(nOld > 0, 'got empty signal');
    realFlag = isreal(td);
    
    if n == nOld
        return;
    end
    
    fd = fft(td);    
    
    if n > nOld
        n1 = ceil(nOld / 2);
        n2 = nOld - n1;
        fd = [fd(1:n1), zeros(1, n-nOld), fd(end-n2+1:end)];
    else
        n1 = ceil(n/2);
        n2 = n - n1;
        fd = fd([1:n1, end-n2+1:end]);
    end
    td = ifft(fd * n / nOld);
    
    if realFlag
        td = real(td);
    end
end

function td = td2td_delay(td, delay_samples)
    if nargin == 0
        d();
        d('td = td2td_delay(td, delay_samples)');
        d();
        d('    delay a cyclic signal');
        d('td:');
        d('    input signal');
        d('delay_samples: (real)');
        d('    delay to apply on td');
        return;
    end
    if delay_samples == 0 
        return;
    end
    rflag = isreal(td);
    
    n = numel(td);
    
    f = 0:(n - 1);
    f = f + floor(n / 2);
    f = mod(f, n);
    f = f - floor(n / 2);
    
    % delay in terms of signal lengths
    nCyc = delay_samples / n;
    phase = -2 * pi * f * nCyc;
    rot = exp(1i*phase);
    
    td = ifft(fft(td) .* rot);
    if rflag
        td = real(td);
    end
end

% fundamental frequency estimation 
function tPeriod = td_estFundPeriod(s, tMax)
    if nargin == 0
        d();
        d('tPeriod = td_estFundPeriod(s, tMin)');
        d();
        d('    estimates the fundamental period in s');
        d('s:');
        d('    input data');
        d('tMin: (optional)');
        d('    limit estimate to tMax (suppress detection of multiples)');
        d('tPeriod (real)');
        d('    fundamental period estimate');
        return;
    end

    if nargin < 2
        tMax = floor(numel(s)/2);
    end
    
    % === autocorrelation ===
    n = numel(s);
    tmp = fft(s);
    tmp = real(ifft(tmp .* conj(tmp)));
    
    % note: could take logarithm here (cepstrum)
    
    % === periodicity of autocorrelation ===
    tmp = fft(tmp);
    tmp = real(ifft(tmp .* conj(tmp)));
    
    % === blank out dt=0 peak ===
    tmp(1) = 0;
    for ix = 2 : numel(tmp)/2-1
        if tmp(ix) < tmp(ix+1)
            break;
        end
        tmp(ix) = 0;
        tmp(end-ix+2) = 0;
    end
    
    % === pick eligible maxima ===
    xc = tmp;
    tmp = xc(1:tMax);
    xcMax = max(tmp);
    
    % === locate maximum between samples ===
    % f(x):=a+b*x+c*x^2;
    % sol:solve([f(-1)=vm, f(0)=v, f(1)=vp], [a, b, c]);
    % solve(diff(ev(f(x), sol), x, 1)=0, x);
    % diff(ev(f(x), sol), x, 2);
    % Solution: "[x = -(vp-vm)/(2*vp+2*vm-4*v)]"
    % 2nd derivative: vp+vm-2*v
    % Poly fit: 
    % [a = v,b = -(vm-vp)/2,c = (vp+vm-2*v)/2]
    
    ixPeak = find(tmp == xcMax, 1, 'first');

    r = ixPeak-1:ixPeak+1;
    r = mod(r-1, n)+1;
    vm = sqrt(double(xc(r(1))));
    v = sqrt(double(xc(r(2))));
    vp = sqrt(double(xc(r(3))));
    
    x = -(vp-vm)/(2*vp+2*vm-4*v);
    ddx = vp+vm-2*v;
    
    if ((ddx >= 0) || (x < -1) || (x > 1)) 
        disp('warning: no useful poly maximum. Using peak xcorr bin');
        x = 0;
    end
    tPeriod = ixPeak - 1 + x;
end

function s = td2td_dcBlock(s)
    if nargin == 0
        d();
        d('s = td2td_dcBlock(s)');
        d();
        d('    removes average from s');
        d('s:');
        d('    input data');
        return;
    end
    s = s - sum(s) / numel(s);
end

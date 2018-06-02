# libsk61

A library with Octave (Matlab) utility functions

## Spline

Spline fit with optional discontinuities in value and / or slope

See example/sk61example_spline.m

### spline.create(...)

Fit a spline to data

### spline.eval(...)

Evaluate a spline

## td2td

Operations on cyclic signals (td2td for 'timedomain to timedomain')

### td2td.noncausalResampler(...)

Changes sample rate using ideal lowpass filtering

### td2td.delay(...)

delay (cyclic shift) with subsample precision

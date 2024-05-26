rstudioapi::jobRunScript('analysis/hr-doy-regression.R',
                         name = 'HR',
                         workingDir = '.')

rstudioapi::jobRunScript('analysis/speed-doy-regression.R',
                         name = 'speed',
                         workingDir = '.')

rstudioapi::jobRunScript('analysis/diffusion-doy-regression.R',
                         name = 'diffusion',
                         workingDir = '.')

rstudioapi::jobRunScript('analysis/density-doy-regression.R',
                         name = 'density',
                         workingDir = '.')

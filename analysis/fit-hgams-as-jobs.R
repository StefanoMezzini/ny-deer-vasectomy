rstudioapi::jobRunScript('analysis/hr-doy-regression.R',
                         name = 'HR',
                         workingDir = '.')

rstudioapi::jobRunScript('analysis/speed-doy-regression.R',
                         name = 'speed',
                         workingDir = '.')

rstudioapi::jobRunScript('analysis/diffusion-doy-regression.R',
                         name = 'diffusion',
                         workingDir = '.')

rstudioapi::jobRunScript('analysis/excursivity-doy-regression.R',
                         name = 'density',
                         workingDir = '.')

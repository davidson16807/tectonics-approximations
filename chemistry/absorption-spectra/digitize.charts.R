# NOTE: 
# "here" is used to ensure to the correct working directory without using global paths
# "digitize" is used to extract point information from png charts
# install.packages(c('here','digitize'))
# setwd(here::here('research/atmosphere/absorption-spectra'))

#digitize charts manually

# what follows are a bunch of functions for fast, crude, well-behaved approximations

# The Zucconni "bump" function is a convenient way to model bell shaped functions
# it's easy to parameterize, integrate, understand, and compute
bump = function(x, edge0, edge1, height){
    center = (edge1 + edge0) / 2.;
    width = (edge1 - edge0) / 2.;
    offset = (x - center) / width;
    return( height * pmax(1. - offset * offset, 0.) )
}
# "clamp" returns the closest approximation of x that exists between lo and hi
clamp = function(x,lo,hi) {
    return(pmin(pmax(x, lo), hi))
}
# "linearstep" returns a fraction indicating how far in x is between lo and hi
#  it can be used with "mix" to easily interpolate between values
linearstep = function(lo, hi, x){
    return(clamp((x - lo) / (hi - lo), 0.0, 1.0));
}
# "lerp" performs basic Linear piecewise intERPolation:
#  given a list of control points mapping 1d space to 1d scalars, 
#  and a point in 1d space, returns a 1d scalar that maps to the point
lerp = function(control_point_x, control_point_y, x) {
    return(sapply(x,function(xi) {
                outi = control_point_y[1];
                for (i in 2:length(control_point_x)) {
                    outi = mix(outi, control_point_y[i], linearstep(control_point_x[i-1], control_point_x[i], xi));
                }
                return(outi)
            }
        )
    );
}
mix = function(x, y, a) {
    return(x*(1.-a) + y*a);
}


cpp.list = function(text) {
  gsub('+0', '', 
    paste('{ ', 
      paste(text, 'f', sep='', collapse=', '), 
    '},'), 
  fixed=T)
}

cpp.sample.vectors = function(wavenumber, cross.section) {
  cat(
    paste(
      '\n',
      cpp.list(format(wavenumber[order(wavenumber)], scientific=T, digits=3)),
      cpp.list(paste(format(log10(cross.section[order(wavenumber)]), scientific=F, digits=4))),
      '\n',
      sep='\n'
    )
  )
}

r.list = function(text) {
  gsub('+0', '', 
    paste('c( ', 
      paste(text, sep='', collapse=', '), 
    '),'), 
  fixed=T)
}

r.sample.vectors = function(wavenumber, cross.section) {
  cat(
    paste(
      '\n',
      r.list(format(wavenumber[order(wavenumber)], scientific=T, digits=3)),
      r.list(paste(format(log10(cross.section[order(wavenumber)]), scientific=F, digits=4))),
      '\n',
      sep='\n'
    )
  )
}

# assumptions:
# * x axis is in nm on a linear scale
# * y axis is in cm^2/molecule on a log scale
cpp.sample.df = function(data) {
  wavenumber = 1/(data$x * 1e-9)
  cross.section = 1e-4 * 10^(data$y)
  cpp.sample.vectors(wavenumber, cross.section)
}
r.sample.df = function(data) {
  wavenumber = 1/(data$x * 1e-9)
  cross.section = 1e-4 * 10^(data$y)
  r.sample.vectors(wavenumber, cross.section)
}

dalton = 1.66054e-27 # kilograms


x=seq(0, 5e8, length.out = 30000)

#convert x from wave number in 1/cm to wavenumber in 1/meter, and y from log10(cm^2/molecule) to m^2/molecule
charts = lapply( charts,  function(chart) { 
    chart$wavenumber = 100*chart$V1
    chart$cross.section = 1e-4 * chart$V2
    return(chart)
})
# some charts are not provided by HITRAN
# they're there to corroborate HITRAN and fill in gaps
# they require custom logic for unit conversion
charts$`H2O-liquid.png` = digitize::digitize('H2O-liquid.png', x1=-8, x2=-2, y1=-1, y2=9)
charts$`H2O-liquid.png`$wavenumber = 1/(10^(charts$`H2O-liquid.png`$x))
charts$`H2O-liquid.png`$cross.section = 10^(charts$`H2O-liquid.png`$y) / (1000/(18.015*dalton))

charts$`H2O-1.jpg` = digitize::digitize('H2O-1.jpg', x1=0, x2=210, y1=-21, y2=-16)
charts$`H2O-1.jpg`$wavenumber = 1/(charts$`H2O-1.jpg`$x * 1e-9)
charts$`H2O-1.jpg`$cross.section = 1e-4 * 10^(charts$`H2O-1.jpg`$y)

charts$`H2O-2.jpg` = digitize::digitize('H2O-2.jpg', x1=280, x2=430, y1=-26, y2=-23)
charts$`H2O-2.jpg`$wavenumber = 1/(charts$`H2O-2.jpg`$x * 1e-9)
charts$`H2O-2.jpg`$cross.section = 1e-4 * 10^(charts$`H2O-2.jpg`$y)

plot(c(), c(), type='lines', 
    xlim=c(0,3e7), ylim=c(-36,-18))
with(charts$`H2O-liquid.png`, points(wavenumber, log10(cross.section), type='lines', col=1, lty=2))
with(charts$`H2O-1.jpg`, points(wavenumber, log10(cross.section), type='lines', col=2))
with(charts$`H2O-2.jpg`, points(wavenumber, log10(cross.section), type='lines', col=3))

range(charts$`H2O-liquid.png`$wavenumber)
range(charts$`H2O-1.jpg`$wavenumber)
range(charts$`H2O-2.jpg`$wavenumber)

charts$`H2O` = rbind(
  charts$`H2O-1.jpg`, 
  charts$`H2O-liquid.png`[
    (charts$`H2O-liquid.png`$wavenumber < 2380850) | 
    (3081527 < charts$`H2O-liquid.png`$wavenumber & charts$`H2O-liquid.png`$wavenumber < 402170) ,
  ], 
  charts$`H2O-2.jpg`
)
charts$`H2O` = charts$`H2O`[charts$`H2O`$cross.section > 0,]
charts$`H2O` = charts$`H2O`[order(charts$`H2O`$wavenumber),]
charts$`H2O`$cross.section = pmax(charts$`H2O`$cross.section, 1e-35)
lines(charts$`H2O`$wavenumber, log10(charts$`H2O`$cross.section), type='lines', col=4, lty=3)
cpp.sample.vectors(charts$`H2O`$wavenumber, charts$`H2O`$cross.section)


# O3 gas
charts$`O3-1.png` = digitize::digitize('O3-1.png')
charts$`O3-1.png`$wavenumber = 1/(charts$`O3-1.png`$x * 1e-9)
charts$`O3-1.png`$cross.section = 1e-4 * 10^(charts$`O3-1.png`$y)

charts$`O2O3-O2.png` = digitize::digitize('O2O3.png')
charts$`O2O3-O2.png`$wavenumber = 1/(charts$`O2O3-O2.png`$x * 1e-9)
charts$`O2O3-O2.png`$cross.section = 10^(charts$`O2O3-O2.png`$y)

charts$`O2O3-O3.png` = digitize::digitize('O2O3.png')
charts$`O2O3-O3.png`$wavenumber = 1/(charts$`O2O3-O3.png`$x * 1e-9)
charts$`O2O3-O3.png`$cross.section = 10^(charts$`O2O3-O3.png`$y)

charts$`O3-2.png` = digitize::digitize('O3-2.png')
charts$`O3-2.png`$wavenumber = 1/(charts$`O3-2.png`$x * 1e-9)
charts$`O3-2.png`$cross.section = 1e-4 * 10^(charts$`O3-2.png`$y)

with(charts$`O3.txt`, plot(wavenumber, log10(cross.section), main='O3', type='lines', xlim=c(0,1.5e7), ylim=c(-36,-20)))
with(charts$`O2O3-O2.png`, points(wavenumber, log10(cross.section), type='lines', col=3))
with(charts$`O2O3-O3.png`, points(wavenumber, log10(cross.section), type='lines', col=4))
with(charts$`O3-1.png`,    points(wavenumber, log10(cross.section), type='lines', col=5))
with(charts$`O3-2.png`,    points(wavenumber, log10(cross.section), type='lines', col=6))
# O3 estimate
lines(x, lerp(c(0,  2e5,7e5,9e5,1.6e6,2e6,2.5e6,2.8e6,3e6,3.5e6,4.6e6,6e6,7.7e6,1.2e7), 
              c(-28,-26,-31,-28,-24,  -25,-27,  -24.5,-23,-21,  -22.5,-22,-21,  -21  ), x), col=2)


charts$`C2H6-1.jpg` = digitize::digitize('C2H6-1.jpg')
charts$`C2H6-1.jpg`$wavenumber = 1/(charts$`C2H6-1.jpg`$x * 1e-9)
charts$`C2H6-1.jpg`$cross.section = 1e-4 * 10^(charts$`C2H6-1.jpg`$y)

with(charts$`C2H6.txt`,   plot(wavenumber,   log10(cross.section), main='C2H6', type='lines', xlim=c(0,1e9), ylim=c(-36,-20)))
with(charts$`C2H6-1.jpg`, points(wavenumber, log10(cross.section), type='lines', col=3))
lines(x, lerp(c(5.6e6, 7.6e6, 1.2e7, 5.3e7, 1.9e8), 
              c(-35,   -20.6, -20,   -21.5, -22.6), x)
         # * abs(sin(x/(2*170e3))) 
         ,col=2)

charts$`N2O-1.jpg` = digitize::digitize('N2O-1.jpg')
charts$`N2O-1.jpg`$wavenumber = 1/(charts$`N2O-1.jpg`$x * 1e-9)
charts$`N2O-1.jpg`$cross.section = 1e-4 * 10^(charts$`N2O-1.jpg`$y)

charts$`N2O-2.jpg` = digitize::digitize('N2O-2.jpg')
charts$`N2O-2.jpg`$wavenumber = 1/(charts$`N2O-2.jpg`$x * 1e-9)
charts$`N2O-2.jpg`$cross.section = 1e-4 * 10^(charts$`N2O-2.jpg`$y)

with(charts$`N2O.txt`,   plot(wavenumber,   log10(cross.section), main='N2O', type='lines', xlim=c(0,1e9), ylim=c(-36,-20)))
with(charts$`N2O-1.jpg`, points(wavenumber, log10(cross.section), type='lines', col=3))
with(charts$`N2O-2.jpg`, points(wavenumber, log10(cross.section), type='lines', col=3))
lines(x, lerp(c(5e4, 6e4, 2.6e5, 7.7e5, 2.7e6, 7.6e6, 7.8e7, 2.1e8, 4.4e8), 
              c(-35, -25, -24.4, -29,   -35,   -20.4, -21.4, -22.4, -21.8), x)
         # * abs(sin(x/(2*170e3))) 
         ,col=2)

# Our coverage of N2 extends only to UV, but it's not very active in Vis/IR
cpp.sample.df(digitize::digitize('N2-1.jpg', x1=0, x2=130, y1=-24, y2=-14))

# Our coverage of O2 extends only to UV, but it's not very active in Vis/IR
o2.1=digitize::digitize('O2-1.jpg', x1=100, x2=190, y1=-21, y2=-16)
o2.2=digitize::digitize('O2-2.jpg', x1=0, x2=120, y1=-19, y2=-16)
cpp.sample.df(rbind(o2.1,o2.2))

# I am super proud of our CO2 coverage: 
# it extends from 2.5nm to 1mm with only one major gap in Vis 
# (we patch the gap using `co2.gap`, since we know CO2 is not active in Vis)
# Three separate graphs are used to construct CO2, 
# and each has zero overlap so there is little need to worry about jagged lines 
# or inefficient sample point placement.
# Truly only the best for such a common and highly influential gas!
co2.1=digitize::digitize('CO2-1.jpg', x1=0, x2=210, y1=-25, y2=-15)
co2.2=digitize::digitize('CO2-2.png', x1=0.5, x2=4.5, y1=-35, y2=-15)
co2.2$x = co2.2$x * 1000 #convert um to nm
co2.3=digitize::digitize('CO2-2.png', x1=0, x2=2000, y1=-35, y2=-15)
co2.3$x = 1e7/co2.3$x #convert 1/cm to nm
co2.gap = data.frame(x=400, y=-36)
cpp.sample.df(rbind(co2.1,co2.2,co2.3,co2.gap))

# The graph we used for CO2 also includes CH4, 
# so our coverage of CH4 is also extensive, from ~2nm to 1mm.
# Three graphs are used, with virtually zero gaps or overlaps! 
# There is only a small gap from 165nm to 440nm, 
# but absorption there is not appreciable.
ch4.1=digitize::digitize('CH4-1.jpg', x1=400, x2=1100, y1=-28, y2=-22)
ch4.2=digitize::digitize('CH4-2.jpg', x1=0, x2=180, y1=-23, y2=-15)
ch4.3=digitize::digitize('CH4-3.png', x1=0.5, x2=4.5, y1=-35, y2=-15)
ch4.3$x = ch4.2$x * 1000 #convert um to nm
ch4.4=digitize::digitize('CH4-3.png', x1=0, x2=2000, y1=-35, y2=-15)
ch4.4$x = 1e7/ch4.3$x #convert 1/cm to nm
cpp.sample.df(rbind(ch4.1,ch4.2,ch4.3))



# Our coverage of CF4 includes UV and IR, 
# We went out of our way searching for an IR graph since 
# we use CF4 to study pollution and Martian terraformation
cf4.2=digitize::digitize('CF4-2.jpg', x1=1282.8, x2=1283.4, y1=-16, y2=log10(4e-16))
cf4.1=digitize::digitize('CF4-1.jpg', x1=0, x2=130, y1=-21, y2=-16)
cf4.2$x = 1e7/cf4.2$x
cpp.sample.df(rbind(cf4.2, cf4.1))

# what follows are less influential gases, 
# we take what we can get from the UV/Vis spectral atlas
cpp.sample.df(digitize::digitize('NH3-1.jpg', x1=0, x2=240, y1=-22, y2=-16))
n2o.1=digitize::digitize('N2O-1.jpg', x1=0, x2=170, y1=-19, y2=-16)
n2o.2=digitize::digitize('N2O-2.jpg', x1=125, x2=325, y1=-27, y2=-16)
cpp.sample.df(rbind(n2o.1,n2o.2))
cpp.sample.df(digitize::digitize('SO2.jpg', x1=200, x2=400, y1=-25, y2=-15))
cpp.sample.df(digitize::digitize('NO.jpg', x1=0, x2=250, y1=-20, y2=-15))
co.1=digitize::digitize('CO-1.jpg', x1=170, x2=210, y1=-24, y2=-20)
co.2=digitize::digitize('CO-2.jpg', x1=0, x2=180, y1=-22, y2=-16)
cpp.sample.df(rbind(co.1,co.2))
cpp.sample.df(digitize::digitize('C2H6-1.jpg', x1=0, x2=180, y1=-21, y2=-16))
cpp.sample.df(digitize::digitize('CH2O.jpg', x1=220, x2=420, y1=-24, y2=-18))
cpp.sample.df(digitize::digitize('HCN.jpg', x1=60, x2=120, y1=-20, y2=-15))
c6h6.2=digitize::digitize('C6H6-2.jpg', x1=220, x2=280, y1=-23, y2=-17)
c6h6.1=digitize::digitize('C6H6-1.jpg', x1=0, x2=220, y1=-20, y2=-15)
cpp.sample.df(rbind(c6h6.2, c6h6.1))



# save, for the love of god, save
save(charts, file='charts.Rdata')

# NOTE TO SELF: nitrous oxide calibration was messed up during the last digitization effort
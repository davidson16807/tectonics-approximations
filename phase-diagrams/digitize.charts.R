# install.packages(c('digitize'))
# devtools::install_github('jimhester/fstrings')
data = load('data.Rdata')


# returns the center and radius of a circle that circumscribes 3 points
# z1,z2,z3: imaginary numbers indicating 2d sample points
get.circle.of.3.points = function(z1,z2,z3) {
    w = (z3-z1) / (z2-z1)
    c = (z2 - z1)*(w - abs(w)^2)/(2i*Im(w)) + z1
    return(list(
        center=c,
        radius=abs(z1 - c)
    ))
}

digitize.simon.glatzwell.parameters = function(filename, instructions='', Ptransform=identity) {
    print(instructions)
    print('identify points along the vaporization line: ')
    pref = as.numeric(readline(prompt="Enter example distance along x axis: "))
    tref = as.numeric(readline(prompt="Enter example distance along x axis: "))
    data = digitize::digitize(filename)
    data$Tinv = 1/data$x
    data$lnP = log(Ptransform(data$y))
    return(lm(Tinv~lnP, data=data, na.action=na.omit)$coefficients)
}
digitize.clapeyron.parameters = function(filename, instructions='', Ptransform=identity) {
    print(instructions)
    print('identify points along the vaporization line: ')
    data = digitize::digitize(filename)
    data$Tinv = 1/data$x
    data$lnP = log(Ptransform(data$y))
    return(lm(Tinv~lnP, data=data, na.action=na.omit)$coefficients)
}
digitize.circle = function(filename, instructions='') {
    print(instructions)
    print('identify 3 points along a circle: ')
    data = digitize::digitize(filename)
    dx = as.numeric(readline(prompt="Enter example distance along x axis: "))
    dy = as.numeric(readline(prompt="Enter distance along the y axis that's the same size on the screen: "))
    circle = get.circle.of.3.points(
        data[1,]$x/dx + data[1,]$y/dy*1i,
        data[2,]$x/dx + data[2,]$y/dy*1i,
        data[3,]$x/dx + data[3,]$y/dy*1i
    )
    return(list(
        center = circle$center,
        radius = circle$radius,
        scalex = dx,
        scaley = dy
    ))
}

digitize.engineering.toolbox = function(filename, instructions='') {
    print(instructions)
    print('identify 3 points along a circle: ')
    data = digitize::digitize(filename)
    dx = as.numeric(readline(prompt="Enter example distance along x axis: "))
    dy = as.numeric(readline(prompt="Enter distance along the y axis that's the same size on the screen: "))
    circle = get.circle.of.3.points(
        data[1,]$x/dx + data[1,]$y/dy*1i,
        data[2,]$x/dx + data[2,]$y/dy*1i,
        data[3,]$x/dx + data[3,]$y/dy*1i
    )
    diagram = list()
    diagram$solidus = data[4:nrow(data),]
    diagram$liquidus = list(
        center = circle$center,
        radius = circle$radius,
        scalex = dx,
        scaley = dy
    )
    return(diagram)
}

# EARTH-LIKE MINERALS
data=list()
data$facies = list()
data$facies$zeolite              = digitize.circle('facies.png')
data$facies$blueschist           = digitize.circle('facies.png')
data$facies$ecologite            = digitize.circle('facies.png')
data$facies$hornfels             = digitize.circle('facies.png')
data$facies$prehnite_pumpellyte  = digitize.circle('facies.png')
data$facies$greenschist          = digitize.circle('facies.png')
data$facies$epidote_amphibiolite = digitize.circle('facies.png')
data$facies$amphibiolite         = digitize.circle('facies.png')

# MOST COMMON PURE ELEMENTS, ORDERED BY ABUNDANCE
data$H = digitize.engineering.toolbox('H.jpg')
data$He = digitize.engineering.toolbox('He.jpeg')
data$N = digitize.engineering.toolbox('N.jpg')

data$C = list()
data$C$vaporus  = digitize::digitize('C-2.png')
data$C$liquidus = digitize::digitize('C-2.png')
data$C$solidus  = digitize::digitize('C-2.png')

data$O = digitize.engineering.toolbox('O-2.png')
data$Ar = digitize.engineering.toolbox('Ar-2.png')

# NOTE: x and y are reversed in this diagram
data$Fe = list()
data$Fe$liquidus1 = digitize.circle('Fe-2.png')
data$Fe$liquidus2 = digitize::digitize('Fe-2.png')
data$Fe$delta = digitize::digitize('Fe-2.png')
data$Fe$gamma = digitize::digitize('Fe-2.png')

# COMMON HYDROGEN COMPOUNDS, ORDERED BY ABUNDANCE OF BOUNDED ELEMENT
data$CH4 = digitize.engineering.toolbox('CH4.jpg')
data$NH3 = digitize.engineering.toolbox('NH3.jpg')
data$C2H6 = digitize.engineering.toolbox('C2H6.png')

data$H2O = list()
data$H2O$liquidus = digitize.circle('H2O-4.png')
data$H2O$solidus  = digitize::digitize('H2O-4.png')

# COMMON OXYGEN COMPOUNDS, ORDERED BY ABUNDANCE OF BOUNDED ELEMENT
data$CO2 = digitize.engineering.toolbox('CO2.jpg')

data$CaCO3 = list()
data$CaCO3$liquidus = digitize.circle('CaCO3.png')
data$CaCO3$solidus  = digitize::digitize('CaCO3.png')

data$H$clapeyron = digitize.clapeyron.parameters('H.jpg')
data$He$clapeyron = digitize.clapeyron.parameters('He.jpeg', Ptransform=function(x){10^x})
data$C$clapeyron = digitize.clapeyron.parameters('C-2.png', Ptransform=function(x){10^x})
data$C$clapeyron2 = digitize.clapeyron.parameters('C-2.png', Ptransform=function(x){10^x})
data$N$clapeyron = digitize.clapeyron.parameters('N.jpg')
data$O$clapeyron = digitize.clapeyron.parameters('O-2.png')
data$Ar$clapeyron = digitize.clapeyron.parameters('Ar-2.png')

data$Fe$clapeyron = digitize::digitize('Fe-2.png')
data$Fe$clapeyron$Tinv = 1/data$Fe$clapeyron$y
data$Fe$clapeyron$lnP = log(10^data$Fe$clapeyron$x)
data$Fe$clapeyron = lm(Tinv~lnP, data=data$Fe$clapeyron, na.action=na.omit)$coefficients

data$CH4$clapeyron = digitize.clapeyron.parameters('CH4.jpg')
data$NH3$clapeyron = digitize.clapeyron.parameters('NH3.jpg')
data$C2H6$clapeyron = digitize.clapeyron.parameters('C2H6.png')
data$H2O$clapeyron = digitize.clapeyron.parameters('H2O-4.png', Ptransform=function(x) {10^x})
data$H2O$clapeyron2 = digitize.clapeyron.parameters('H2O-4.png', Ptransform=function(x) {10^x})
data$CO2$clapeyron = digitize.clapeyron.parameters('CO2.jpg')
data$CO$clapeyron = digitize.clapeyron.parameters('CO.png')
data$SiO2 = digitize::digitize('SiO2-3.png')

save(data, file='facies.Rdata')


format.C = function(data) {
    data$t = data$x
    data$p = 10^data$y
    return(data)
}
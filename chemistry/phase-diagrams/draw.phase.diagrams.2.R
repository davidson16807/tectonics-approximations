



    get_basic_phase = function(
        p, # current pressure
        t, # current temperature
        p0,# triple point pressure
        t0,# triple point temperature
        pc,# critical point pressure
        tc,# critical point temperature
        pf,# freezing point pressure
        tf,# freezing point temperature
        L, # latent heat of vaporization at boiling point
        M, # molar mass
        a, # simon glatzel slope
        c, # simon glatzel exponent
        has_simon_glatzel_parameters,
        has_critical_point
    ){
        Ru = 8.3144598
        R = Ru/M
        m = -(R/L)
        b = if(has_critical_point) (1/tc - m * log(pc)) else (1/t0 - m * log(p0))
        A = (tf-t0) / (pf-p0)
        if (t > tc && p > pc && has_critical_point)
        {
            return(3);
        }
        else if ( t > 1/(m*log(p)+b) )
        {
            return(2);
        }
        else if ( t > ( t0 * max((p-p0)/a + 1, 0)^(1/c) ) && has_simon_glatzel_parameters )
        {
            return(1);
        }
        else if ( t > t0 + A*p && !has_simon_glatzel_parameters ) # linear relation, for when simon glatzel exponent can't be determined
        {
            return(1);
        }
        else
        {
            return(0);
        }
    }
p=rep(seq(0.0,3e6,length=50), times=50)
t=rep(seq(0.0,4e3,length=50), each=50)

                                                 # p,    t,    p0,       t0,       pc,       tc,      pf,     tf,      L (J/kg)   M,          a,        c
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 7.04e3  , 13.84   , 1300e3  , 33.20  , 101e3 , 13.99  , 447000   , 0.002016  , 274.22e5, 1.74407 , TRUE  ,  TRUE  )}) # Hydrogen          
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 5.048e3 , 2.1768  , 227e3   , 5.19   , 101e3 , 0.95   , 20700    , 0.004     , 50.96e5 , 1.5602  , TRUE  ,  TRUE  )}) # Helium            
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 12.6e3  , 63.1604 , 3390e3  , 126.2  , 101e3 , 63.15  , 199000   , 0.028     , 1607e5  , 1.7910  , TRUE  ,  TRUE  )}) # Nitrogen          
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 0.152e3 , 54.383  , 5050e3  , 154.6  , 101e3 , 54.36  , 213000   , 0.032     , 2733e5  , 1.7425  , TRUE  ,  TRUE  )}) # Oxygen            
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 68.9e3  , 83.812  , 4870e3  , 150.8  , 101e3 , 83.8   , 163000   , 0.040     , 2114e5  , 1.593   , TRUE  ,  TRUE  )}) # Argon 
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 1       , 1811    , 5100e5  , 8500   , 101e3 , 1811   , 6340000  , 0.056     , 1.07e11 , 1.76    , TRUE  ,  TRUE  )}) # Iron                 lower bound for Pc is based on critical point of gold
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 11.7e3  , 90.66   , 4640e3  , 190.8  , 101e3 , 90.68  , 511000   , 0.016     , 2080e5  , 1.698   , TRUE  ,  TRUE  )}) # Methane           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 6.060e3 , 195.30  , 11280e3 , 405.5  , 101e3 , 195.4  , 1369000  , 0.017     , 5270e5  , 4.3     , TRUE  ,  TRUE  )}) # Ammonia           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 0.0014e3, 90.35   , 4890e3  , 305.3  , 101e3 , 90.36  , 487645.2 , 0.03      , 0       , 0       , FALSE ,  TRUE  )}) # Ethane            
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 15.37e3 , 68.146  , 3500e3  , 132.85 , 101e3 , 68.15  , 216000   , 0.028     , 0       , 0       , FALSE ,  TRUE  )}) # Carbon Monoxide   
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 517e3   , 216.56  , 7380e3  , 304.19 , 101e3 , 216.6  , 205000   , 0.044     , 4000e5  , 2.60    , TRUE  ,  TRUE  )}) # Carbon Dioxide    
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 87.85e3 , 182.34  , 7240e3  , 309.5  , 101e3 , 182.2  , 376141.6 , 0.044     , 0       , 0       , FALSE ,  TRUE  )}) # Nitric Oxide     
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 1.67e3  , 197.69  , 7873e3  , 430.65 , 101e3 , 200.2  , 381464   , 0.064     , 0       , 0       , FALSE ,  TRUE  )}) # Sulfur dioxide    
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 0.1012e3, 89.54   , 3793e3  , 227.6  , 101e3 , 89.15  , 137000   , 0.088     , 0       , 0       , FALSE ,  TRUE  )}) # Tetraflouromethane
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 0.153e3 , 259.86  , 5400e3  , 456.65 , 101e3 , 259.7  , 1030000  , 0.026019  , 3080e5  , 3.6     , TRUE  ,  TRUE  )}) # Cyanide           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 4.3e-4  , 150     , 6300e3  , 514    , 101e3 , 159.2  , 855000   , 0.046     , 10600e5 , 1.61    , TRUE  ,  TRUE  )}) # Ethanol           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 71549032, 155.10  , 6788e3  , 410.3  , 101e3 , 181    , 393296   , 0.030     , 0       , 0       , FALSE ,  TRUE  )}) # Formaldehyde      
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 2.2e3   , 281.4   , 5810e3  , 588    , 101e3 , 281.5  , 437000   , 0.046     , 4100e5  , 5.2     , TRUE  ,  TRUE  )}) # Formic Acid       
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], 0.0003  , 1983    , 1.7e8   , 5300   , 101e3 , 1983   , 11770e3  , 0.06008   , 0       , 0       , FALSE ,  TRUE  )}) # Silica 
#phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],         ,         ,         ,        ,       ,        , 8205e3   , 0.040304  ,         ,         , FALSE ,  TRUE  )}) # Magnesium Oxide
    plot(t,p,col=phase+1, pch=15)
                                                                                                                                                                                                         



# carbon has two solid phases
# one of them (diamond) is economically relevant
# so we give special attention to it

    get_carbon_phase = function(
        p,  # current pressure
        t   # current temperature
    ){
        Ru = (8.3144598);           # universal gas constant
        M  = (0.012);               # molar mass, kg/mol
        Lv = (29650000);            # specific latent heat of vaporization (J/kg)
        p0 = (0.2e6);               # triple point pressure (Pa)
        t0 = (4600);                # triple point temperature (K)
        R  = (Ru/M);                # individual gas constant
        mv = (-R/Lv);               # slope of clapeyron equation for vaporus
        bv = (1/t0 - mv * log(p0)); # intercept for clapeyron equation for vaporus
        ml = (-3.411e-5);           # slope of clapeyron equation repurposed for liquidus, estimated from phase diagram
        bl = (1e-3);                # intercept for clapeyron equation repurposed for liquidus, estimted from phase diagram
        if ( t > 1/(mv*log(p)+bv) && p<1e9)
        {
            return(2);
        }
        else if ( t < t0 && p < 1e10 )
        {
            return(0);
        }
        else if ( t > 1/(ml*log(p)+bl))
        {
            return(1);
        }
        else 
        {
            return(3);
        }
    }


p=rep(10^(seq(0.0,11,length=50)), times=50)
t=rep(seq(0.0,1e4,length=50), each=50)
phase=sapply(seq(p), function(i){ get_carbon_phase( p[i], t[i] )}) 
plot(t,p,col=phase+1, pch=15, log='y')



# Water occurs often in the universe, its one of the simplest compounds that can be made 
# from two of the most abundant elements in the universe
# It has many kinds of solid phases. We're used to ice1h, but ice1c occurs in clouds on earth,
# and ice 7 is common on other planets and moons, so we give its phase diagram special attention.
# Since it takes a long time to go back and revisit phase diagrams if an alteration is needed,
# We will attempt to represent all phases that are currently known to science.
# We don't expect we'll be using these, we just want to make sure we don't have to go through this again.
# Nevertheless, some phases are not very common so they will be represented only in a crude way. 
get_water_phase = function(
    p,  # current pressure
    t   # current temperature
){
    Ru = (8.3144598);          # universal gas constant
    M  = (0.0180153);          # molar mass, kg/mol
    L  = (22.6e5);             # specific latent heat of vaporization (J/kg)
    p0 = (0.6116e3);           # triple point pressure (Pa)
    t0 = (273.15);             # triple point temperature (K)
    R  = (Ru/M);               # individual gas constant
    m  = (-R/L);               # slope of clapeyron equation for vaporus
    b  = (1/t0 - m * log(p0)); # intercept for clapeyron equation for vaporus
    tc = (647.096);            # critical point temperature (K) 
    pc = (22.064e6);           # critical point pressure (Pa)
    a3 = (7070e5);             # slope for the simon and glatzwell equation
    c3 = (4.46);               # intercept for the simon and glatzwell equation
    if ( t > 1/(m*log(p)+b) && p<1e9)
    {
        if (t > tc && p > pc)
        {
            return(3); # supercritical
        }
        else
        {
            return(2); # liquid
        }
    }
    else if ( t > ( t0 * max((p-p0)/a3 + 1, 0)^(1/c3) ) )
    {
        if (t > tc && p > pc)
        {
            return(3); # supercritical
        }
        else
        {
            return(1); # liquid
        }
    }
    else if (p < 209.9e6 && t < 73.15)
    {
        return(2); #11
    }
    else if (p < 209.9e6 && t < 173.15)
    {
        return(3); #1c
    }
    else if (p < 209.9e6)
    {
        return(0); #1h
    }
    else if (p < 632.4e6 && t < 170)
    {
        return(9); #9
    }
    else if (p < 350.1e6 && t > 238.5)
    {
        return(3); #3
    }
    else if (350.1e6 < p && p < 632.4e6 && t > 218)
    {
        return(5); #5
    }
    else if (p < 632.4e6)
    {
        return(2); #2
    }
    else if (p < 2.216e9 && t < 130)
    {
        return(15); #15
    }
    else if (p < 2.216e9)
    {
        return(6); #6
    }
    else if (p < 62e9 && (t-100)/(278-100) < (p-62e9)/(2.1e9-62e9))
    {
        return(1); #8
    }
    else if (p < 62e9)
    {
        return(7); #7
    }
    else if (p < 350e9)
    {
        return(10); #10
    }
    else 
    {
        return(11);
    }
}

p=rep(10^(seq(0,12,length=75)), times=75)
t=rep(seq(0.0,800,length=75), each=75)
phase=sapply(seq(p), function(i){ get_water_phase( p[i], t[i] )})
plot(t,p,col=phase+1, pch=15, log='y')






    # Metamorphic facies are like phases of matter for rocks.
    # Rocks are complex assortments of mineral grains. 
    # Each mineral has its own phase diagram, 
    # but there are many kinds of minerals and their phase diagrams are poorly documented. 
    # However the phase diagrams for minerals share some consistencies,
    # and these consistencies allow us to make predictions about what an earth-like rock will look like
    # when submitted to a given pressure and temperature.
    # This will apply to any rock that is composed of minerals commonly found on earth.
    get_metamorphic_facies = function(p, t)
    {
        pt = c(p,t) / c(2e8, 100);
        if ( sqrt(sum((pt_scaled-c(0,2.65))^2)) < 0.011 )
        {
            return(0) # igneous_or_sediment;
        }
        else if ( sqrt(sum((pt_scaled-c(0,2.65))^2)) < 2.0 )
        {
            return(1) # sedimentary;
        }
        else if( sqrt(sum((pt_scaled-c(0,2.65))^2)) < 3.0 )
        {
            return(2) # zeolite;
        }
        else if( sqrt(sum((pt_scaled-c(9.1,-4.8))^2)) < 11.5 )
        {
            return(4) # blueschist;
        }
        else if( sqrt(sum((pt_scaled-c(13.5,7.1))^2)) < 7.4 )
        {
            return(5) # eclogite;
        }
        else if( sqrt(sum((pt_scaled-c(-9.2,7.6))^2)) < 10.5 )
        {
            return(6) # hornfels;
        }
        else if( sqrt(sum((pt_scaled-c(0.9,-9.2))^2)) < 15.4 )
        {
            return(7) # prehnite_pumpellyte;
        }
        else if( sqrt(sum((pt_scaled-c(5.9,-8.15))^2)) < 15.6 )
        {
            return(2) # greenschist;
        }
        else if( sqrt(sum((pt_scaled-c(6.6, -11.8))^2)) < 19.8 )
        {
            return(9) # epidote_amphibiolite;
        }
        else if( sqrt(sum((pt_scaled-c(6.1,-0.75))^2)) < 10.5 )
        {
            return(10) # amphibolite;
        }
        else
        {
            return(11) # granulite;
        }
    }


p=rep((seq(0,20e8,length=75)), times=75)
t=rep(seq(0,1000,length=75), each=75)
phase=sapply(seq(p), function(i){ get_metamorphic_facies( p[i], t[i]+273.15 )})
plot(t,-p,col=phase+1, pch=15, xlab='temperature (C)')


    get_clapeyron_slope = function(L,M) {
        R_u = 8.3144598
        R = R_u/M
        m = -(R/L)
        return(m)
    }

    get_clapeyron_intercept = function(m, p, t) {
        b = (1/t) - m*log(p)
        return(b)
    }

    get_basic_phase = function(
        p, # current pressure
        t, # current temperature
        m, # clapeyron slope
        b, # clapeyron intercept
        a, # simon glatzel slope
        c, # simon glatzel exponent
        p0,# triple point pressure
        t0,# triple point temperature
        pc,# critical point pressure
        tc,# critical point temperature
        has_simon_glatzel_exponent,
        has_critical_point
    ){
        if (t > tc && p > pc && has_critical_point)
        {
            return(3);
        }
        else if ( t > 1/(m*log(p)+b) )
        {
            return(2);
        }
        else if ( p < p0 )
        {
            return(0);
        }
        else if ( t > ( t0 * ((p-p0)/a + 1)^(1/c) ) && has_simon_glatzel_exponent )
        {
            return(1);
        }
        else if ( t > t0 + (p-p0)/a ) # linear relation, for when simon glatzel exponent can't be determined
        {
            return(1);
        }
        else
        {
            return(0);
        }
    }


p=rep(seq(0,5e6,length=50), times=50)
t=rep(seq(0,400,length=50), each=50)

# Substance           T0      P0(kPa)    a(bars)  c        L(J/kg)   M(kg/mol) m = -R/L  b       Tc     Pc(kPa)   
# Hydrogen            13.84   7.04       274.22   1.74407  455000    0.002016  -0.009064 0.1564  33.20  1300    
#                                                  p     t      m           b        a         c        p0         t0        pc          tc  
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.009064  , 0.1564 , 274.22  , 1.7441 , 7.04e3   , 13.84   , 1300e3    , 33.20  , TRUE , TRUE ) }) #Helium                                           
plot(t,p,col=phase+1, pch=15)
phase=sapply(seq(p), function(i){ g.c_basic_phase( p[i], t[i], -0.001379  , 0.0287 , 1606.5e5, 1.7910 , 12.6e3   , 63.1604 , 3390e3    , 126.2  , TRUE , TRUE ) }) #Nitrogen                                         
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.080453  , 1.1557 , 50.96e5 , 1.5602 , 5.048e3  , 2.1768  , 227e3     , 5.19   , TRUE , TRUE ) }) #Helium                                           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.002121  , 0.0384 , 2732.9e5, 1.7425 , 0.152e3  , 54.383  , 5050e3    , 154.6  , TRUE , TRUE ) }) #Oxygen                                           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.001311  , 0.0268 , 2114e5  , 1.593  , 68.9e3   , 83.812  , 4870e3    , 150.8  , TRUE , TRUE ) }) #Argon                                            
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -2.81e-5   , 5.2e-4 , 1.07e11 , 1.76   , 0         , 1805    , 510000e3 , 8500   , FALSE, TRUE ) }) #Iron              lower bound for Pc is based on critical point of gold                                 
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000140  , 0.0085 , 36.629e5, 0      , 11.7e3   , 90.66   , 4640e3    , 190.8  , FALSE, TRUE ) }) #Methane                                         
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000386  , 0.0087 , 5270e5  , 4.3    , 6.060e3  , 195.30  , 11280e3   , 405.5  , TRUE , TRUE ) }) #Ammonia                                          
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000190  , 0.0048 , -3952e5 , 9.0    , 0.6116e3 , 273.15  , 22060e3   , 647.1  , TRUE , TRUE ) }) #Water I(vapor) 
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.001277  , 0.0223 , 61.4e5, , 0      , 8e-1     , 90.35   , 4890e3    , 305.3  , FALSE, TRUE ) }) #Ethane                                         
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.001375  , 0.0282 , 45.27e5 , 0      , 15.37e3  , 68.146  , 3.49e6    , 132.85 , FALSE, TRUE ) }) #Carbon Monoxide   m and b are calculated from latent heat and may be off
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000465  , 0.0107 , 4000e5  , 2.60   , 517e3    , 216.56  , 7380e3    , 304.19 , TRUE , TRUE ) }) #Carbon Dioxide                                  
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000601  , 0.0031 , 82.29e5 , 0      , 21.92e3  , 109.50  , 6486e3    , 180.15 , FALSE, TRUE ) }) #Nitric Oxide
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000341  , 0.0025 , 0       , 0      , 1.67e3   , 197.69  , 7873e3    , 430.65 , TRUE , TRUE ) }) #Sulfur dioxide
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -0.000689  , 0.0080 , 0       , 0      , 26e3     , 173.08  , 4059e3    , 374.21 , TRUE , TRUE ) }) #Tetraflouromethane
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 0       , 0      ,          , 259.86  , 5400e3    , 456.65 , TRUE , TRUE ) }) #Cyanide           
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 0       , 0      ,          , 155.10  , 6788e3    , 137.2  , TRUE , TRUE ) }) #Formaldehyde                                  
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 0       , 0      , 4.3e-4   , 150     , 6300e3    , 514    , TRUE , TRUE ) }) #Ethanol                                       
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 0       , 0      , 2.2e3    , 281.4   ,           ,        , TRUE , FALSE) }) #Formic Acid

phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -2.61e-5   , 6.3e-4 ,         ,        , 10132e3  , 4765    ,           ,        , TRUE , FALSE) }) #Carbon(vapor)                                      
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i], -e.48e-5   , 1.0e-3 ,         ,        , 10132e3  , 4765    ,           ,        , TRUE , FALSE) }) #Carbon(liquid)                                    
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 620e5   , 60     ,          , 251.20  ,           ,        , TRUE , FALSE) }) #Water III                                       
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 4100e5  , 8.1    ,          , 256.2   ,           ,        , TRUE , FALSE) }) #Water V                                         
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 7070e5  , 4.46   ,          , 273.32  ,           ,        , TRUE , FALSE) }) #Water VI                                        
phase=sapply(seq(p), function(i){ get_basic_phase( p[i], t[i],            ,        , 12980e5 , 3.11   ,          , 354.80  ,           ,        , TRUE , FALSE) }) #Water VII                                        


get_clapeyron_intercept(get_clapeyron_slope(216000, 0.028), 15.37, 68.146)
get_clapeyron_slope(216000, 0.028)
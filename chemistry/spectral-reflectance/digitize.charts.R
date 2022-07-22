
calcite.lo = digitize::digitize('minerals/minerals05')
calcite.hi = digitize::digitize('minerals/calcite-dolomite01')
calcite = digitize::digitize('minerals/minerals28')
aragonite = digitize::digitize('minerals/minerals28')

olivine = digitize::digitize('minerals/olivine-plagioclase-pyroxene-ice.png')
plagioclase = digitize::digitize('minerals/olivine-plagioclase-pyroxene-ice.png')
pyroxene = digitize::digitize('minerals/olivine-plagioclase-pyroxene-ice.png')

pyroxene.hi = digitize::digitize('minerals/pyroxene01')

quartz = digitize::digitize('minerals/quartz05')
cristobalite = digitize::digitize('minerals/quartz05')

hematite = digitize::digitize('minerals/hematite01')
goethite = digitize::digitize('minerals/hematite-goethite')


nacl = digitize::digitize('minerals/salt01')



h2o.liquid = digitize::digitize('cryominerals/ice09')
h2o.ice digitize::digitize('cryominerals/ice-water02')

ch4.liquid.lo = digitize::digitize('cryominerals/methane01')
ch4.liquid.hi = digitize::digitize('cryominerals/methane01')
ch4.ice.lo = digitize::digitize('cryominerals/methane01')
ch4.ice.hi = digitize::digitize('cryominerals/methane01')

n2.ice.alpha = digitize::digitize('cryominerals/n2-01')
n2.ice.beta = digitize::digitize('cryominerals/n2-01')

co2.ice.hi = digitize::digitize('cryominerals/co2-06')

nh3.ice = digitize::digitize('cryominerals/nh3-02')

ethane.ice = digitize::digitize('cryominerals/ethane-01')
ethane.liquid = digitize::digitize('cryominerals/ethane-01')



h2o.ice.alpha = digitize::digitize('cryominerals/nh3-01.png')
h2o.ice.k = digitize::digitize('cryominerals/cryocrystals.png')

co2.ice.alpha = digitize::digitize('cryominerals/nh3-01.png')
co2.ice.k = digitize::digitize('cryominerals/cryocrystals.png')

ch4.ice.alpha = digitize::digitize('cryominerals/nh3-01.png')
ch4.ice.k = digitize::digitize('cryominerals/cryocrystals.png')

so2.ice.alpha = digitize::digitize('cryominerals/nh3-01.png')
so2.ice.k = digitize::digitize('cryominerals/cryocrystals.png')

nh3.ice.k = digitize::digitize('cryominerals/cryocrystals.png')

co.ice.alpha = digitize::digitize('cryominerals/nh3-01.png')


# we've already saved results from the above code into csv, 
# so you can salphaip the process above by calling the code below
# x is wavenumber, y is log10 absorption coefficient in 1/cm
h2o.ice.k = read.csv('h2o.ice.k.csv')
h2o.ice.alpha = read.csv('h2o.ice.alpha.csv')
co2.ice.k = read.csv('co2.ice.k.csv')
co2.ice.alpha = read.csv('co2.ice.alpha.csv')
ch4.ice.k = read.csv('ch4.ice.k.csv')
ch4.ice.alpha = read.csv('ch4.ice.alpha.csv')
so2.ice.k = read.csv('so2.ice.k.csv')
so2.ice.alpha = read.csv('so2.ice.alpha.csv')
nh3.ice.k = read.csv('nh3.ice.k.csv')
# nh3.ice.alpha = read.csv('nh3.ice.alpha.csv')
# co.ice.k = read.csv('co.ice.k.csv')
co.ice.alpha = read.csv('co.ice.alpha.csv')



benzene = digitize::digitize('cryominerals/hydrocarbons.png')


print.vectors = function(x, y) {
	print('                get_absorption_coefficient_function_from_reflectance_at_wavelengths')
    print('                    (si::nanometer, 100.0 * si::micrometer,')
	print(paste('                    std::vector<double>{', paste(x, collapse=', '), '},'))
	print(paste('                    std::vector<double>{', paste(y, collapse=', '), '}),'))
}


print.vectors(format(co2.ice.alpha$x, digits=5, scientific=F, width=7), format(co2.ice.alpha$y, digits=3, scientific=F, width=7))
print.vectors(format(ch4.ice.alpha$x, digits=5, width=7), format(ch4.ice.alpha$y, digits=2, width=7))
print.vectors(format(so2.ice.alpha$x, digits=5, width=7), format(so2.ice.alpha$y, digits=2, width=7))
print.vectors(format(co.ice.alpha$x, digits=5, width=7), format(co.ice.alpha$y, digits=2, width=7))
print.vectors(format(h2o.ice.k$x, digits=5, width=7), format(h2o.ice.k$y, digits=2, width=7))
print.vectors(format(co2.ice.k$x, digits=5, width=7), format(co2.ice.k$y, digits=3, width=7))
print.vectors(format(ch4.ice.k$x, digits=5, width=7), format(ch4.ice.k$y, digits=3, width=7))
print.vectors(format(so2.ice.k$x, digits=5, width=7), format(so2.ice.k$y, digits=3, width=7))
print.vectors(format(nh3.ice.k$x, digits=5, width=7), format(nh3.ice.k$y, digits=3, width=7))


wavenumbers = function(compound, diameter, x, y) {
	print(paste('                get_absorption_coefficient_function_from_reflectance_at_wavenumbers // for', compound))
    print(paste('                    (1.0 / si::centimeter, ', diameter, '*si::micrometer, '))
	print(paste('                    std::vector<double>{', paste(x, collapse=', '), '},'))
	print(paste('                    std::vector<double>{', paste(y, collapse=', '), '}),'))
}
wavelengths = function(compound, diameter, x, y) {
	print(paste('                get_absorption_coefficient_function_from_reflectance_at_wavelengths // for', compound))
    print(paste('                    (si::micrometer, ', diameter, '*si::micrometer, '))
	print(paste('                    std::vector<double>{', paste(x, collapse=', '), '},'))
	print(paste('                    std::vector<double>{', paste(y, collapse=', '), '}),'))
}




# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.1B_ASDFRa_AREF.gif.png')  
# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.3B_NIC4aau_RREF.gif.png')
# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.2B_ASDFRa_AREF.gif.png')  
# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.3B_NIC4aau_RREF_wavenumber.gif.png')
# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.3B_ASDFRa_AREF.gif.png')  
# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.4B_ASDFRa_AREF.gif.png')
# halite.reflectance = digitize::digitize('usgs/splib07a_Halite_HS433.6_ASDFRa_AREF.gif.png')
halite.reflectance.lo = digitize::digitize('usgs/splib07a_Halite_HS433.3B_BECKa_AREF.gif.png')
halite.reflectance.hi = digitize::digitize('usgs/splib07a_Halite_HS433.3B_NIC4aau_RREF.gif.png')
halite.reflectance = rbind(halite.reflectance.lo, halite.reflectance.hi)
write.csv(halite.reflectance, 'usgs/halite.reflectance.csv')
wavelengths('halite', mean(c(74,250)), format(halite.reflectance$x, digits=5, width=7), format(halite.reflectance$y, digits=3, width=7) )

# corundum.reflectance = digitize::digitize('usgs/splib07a_Corundum_HS283.1B_ASDFRc_AREF.gif.png')  
# corundum.reflectance = digitize::digitize('usgs/splib07a_Corundum_HS283.2B_ASDFRc_AREF.gif.png')  
# corundum.reflectance = digitize::digitize('usgs/splib07a_Corundum_HS283.4B_ASDFRc_AREF.gif.png')
corundum.reflectance = digitize::digitize('usgs/splib07a_Corundum_HS283.3B_BECKc_AREF.gif.png')
write.csv(corundum.reflectance, 'usgs/corundum.reflectance.csv')
wavelengths('corundum', mean(c(74,250)), format(corundum.reflectance$x, digits=5, width=7), format(corundum.reflectance$y, digits=3, width=7) )

# apatite.reflectance = digitize::digitize('splib07a_Chlorapatite_REE_WS422_xtl_ASDFRa_AREF.gif.png')   
# apatite.reflectance = digitize::digitize('splib07a_Chlorapatite_WS423_NIC4bb_RREF.gif.png') 
# apatite.reflectance = digitize::digitize('splib07a_Chlorapatite_WS423_ASDFRc_AREF.gif.png')           
# apatite.reflectance = digitize::digitize('splib07a_Chlorapatite_WS423_NIC4bb_RREF_wavenumber.gif.png') 
# apatite.reflectance = digitize::digitize('splib07a_Chlorapatite_WS423_BECKc_AREF.gif.png') 
# apatite.reflectance = digitize::digitize('usgs/splib07a_Hydroxyl-Apatite_WS425_BECKb_AREF.gif.png')   
# apatite.reflectance = digitize::digitize('usgs/splib07a_Hydroxyl-Apatite_WS425_NIC4bb_RREF_wavenumber.gif.png')
# apatite.reflectance = digitize::digitize('usgs/splib07a_Hydroxyl-Apatite_WS425_NIC4bb_RREF.gif.png')
# apatite.reflectance.lo = digitize::digitize('usgs/splib07a_Hydroxyl-Apatite_WS425_BECKb_AREF.gif.png')
# apatite.reflectance.hi = digitize::digitize('usgs/splib07a_Hydroxyl-Apatite_WS425_NIC4bb_RREF.gif.png')
apatite.reflectance = rbind(apatite.reflectance.lo, apatite.reflectance.hi)
write.csv(apatite.reflectance, 'usgs/apatite.reflectance.csv')
# # grain size is not known
# wavelengths('apatite', ), format(apatite.reflectance$x, digits=5, width=7), format(apatite.reflectance$y, digits=3, width=7) )

# carbon.reflectance = digitize::digitize('usgs/splib07a_Carbon_Black_GDS68_ASDFRa_AREF.gif.png')  
carbon.reflectance = digitize::digitize('usgs/splib07a_Carbon_Black_GDS68_sm.ap._BECKa_AREF.gif.png')
write.csv(carbon.reflectance, 'usgs/carbon.reflectance.csv')
wavelengths('carbon', 0.17, format(carbon.reflectance$x, digits=5, width=7), format(carbon.reflectance$y, digits=3, width=7) )

# calcite.reflectance = digitize::digitize('usgs/splib07a_Calcite_CO2004_BECKb_AREF.gif.png')                 
# calcite.reflectance = digitize::digitize('usgs/splib07a_Calcite_WS272_ASDNGa_AREF.gif.png')  
# calcite.reflectance = digitize::digitize('usgs/splib07a_Calcite_GDS304_75-150um_ASDFRb_AREF.gif.png')       
# calcite.reflectance = digitize::digitize('usgs/splib07a_Calcite_HS48.3B_BECKa_AREF.gif.png')                
# calcite.reflectance = digitize::digitize('usgs/splib07a_Calcite_REE-bearing_WS319a_ASDFRb_AREF.gif.png')    
# calcite.reflectance = digitize::digitize('usgs/splib07a_Calcite_WS272_NIC4aaa_RREF_wavenumber.gif.png')  
calcite.reflectance.lo = digitize::digitize('usgs/splib07a_Calcite_WS272_BECKa_AREF.gif.png')  
calcite.reflectance.hi = digitize::digitize('usgs/splib07a_Calcite_WS272_NIC4aaa_RREF.gif.png')  
calcite.reflectance = rbind(calcite.reflectance.lo, calcite.reflectance.hi)
write.csv(calcite.reflectance, 'usgs/calcite.reflectance.csv')
wavelengths('calcite', 410.0, format(calcite.reflectance$x, digits=5, width=7), format(calcite.reflectance$y, digits=3, width=7) )

# quartz.reflectance = digitize::digitize('usgs/splib07a_Quartz_GDS74_Sand_Ottawa_NIC4cb_AREF_wavenumber.gif.png')
# quartz.reflectance = digitize::digitize('usgs/splib07a_Quartz_GDS74_Sand_Ottawa_NIC4cb_AREF_wavenumber.gif.png')
# quartz.reflectance = digitize::digitize('usgs/splib07a_Quartz_HS32.1B_ASDFRc_AREF.gif.png')
# quartz.reflectance = digitize::digitize('usgs/splib07a_Quartz_HS32.2B_ASDFRc_AREF.gif.png')
# quartz.reflectance = digitize::digitize('usgs/splib07a_Quartz_HS32.3B_ASDFRc_AREF.gif.png')
# quartz.reflectance = digitize::digitize('usgs/splib07a_Quartz_HS32.4B_BECKa_AREF.gif.png')
# # We would ideally want GDS74 for both high and low wavelengths, but we are missing grain size for GDS74.
# # We don't have any other samples for Quartz at high wavelengths, so we use one of the other samples 
# # with known grain size for low wavelength and use its size as a representative value for the high wavelength sample
# quartz.reflectance.lo = digitize::digitize('usgs/splib07a_Quartz_GDS74_Sand_Ottawa_BECKc_AREF.gif.png')
quartz.reflectance.lo = digitize::digitize('usgs/splib07a_Quartz_GDS31_0-74um_fr_BECKa_AREF.gif.png')
quartz.reflectance.hi = digitize::digitize('usgs/splib07a_Quartz_GDS74_Sand_Ottawa_NIC4cb_AREF.gif.png')
quartz.reflectance.hi = quartz.reflectance.hi[quartz.reflectance.hi$x > 3,]
quartz.reflectance = rbind(quartz.reflectance.lo, quartz.reflectance.hi)
write.csv(quartz.reflectance, 'usgs/quartz.reflectance.csv')
wavelengths('quartz', mean(c(0, 75)), format(quartz.reflectance$x, digits=5, width=7), format(quartz.reflectance$y, digits=3, width=7) )

# orthoclase.reflectance = digitize::digitize('usgs/splib07a_Orthoclase_NMNH113188_BECKb_AREF.gif.png')     
# orthoclase.reflectance = digitize::digitize('usgs/splib07a_Orthoclase_NMNH142137_Fe_NIC4bb_RREF_wavenumber.gif.png')
orthoclase.reflectance.lo = digitize::digitize('usgs/splib07a_Orthoclase_NMNH142137_Fe_BECKb_AREF.gif.png')  
orthoclase.reflectance.hi = digitize::digitize('usgs/splib07a_Orthoclase_NMNH142137_Fe_NIC4bb_RREF.gif.png')
orthoclase.reflectance.hi = orthoclase.reflectance.hi[orthoclase.reflectance.hi$x > 3,]
orthoclase.reflectance = rbind(orthoclase.reflectance.lo, orthoclase.reflectance.hi)
write.csv(orthoclase.reflectance, 'usgs/orthoclase.reflectance.csv')
wavelengths('orthoclase', mean(c(0, 75)), format(orthoclase.reflectance$x, digits=5, width=7), format(orthoclase.reflectance$y, digits=3, width=7) )

# andesine.reflectance = digitize::digitize('usgs/splib07a_Andesine_HS142.1B_ASDFRc_AREF.gif.png')              
# andesine.reflectance = digitize::digitize('usgs/splib07a_Andesine_HS142.2B_ASDFRc_AREF.gif.png')              
# andesine.reflectance = digitize::digitize('usgs/splib07a_Andesine_HS142.3B_ASDFRc_AREF.gif.png')              
# andesine.reflectance = digitize::digitize('usgs/splib07a_Andesine_HS142.3B_NIC4ccc_RREF_wavenumber.gif.png')  
# andesine.reflectance = digitize::digitize('usgs/splib07a_Andesine_HS142.4B_ASDFRc_AREF.gif.png')
andesine.reflectance.lo = digitize::digitize('usgs/splib07a_Andesine_HS142.3B_BECKc_AREF.gif.png')               
andesine.reflectance.hi = digitize::digitize('usgs/splib07a_Andesine_HS142.3B_NIC4ccc_RREF.gif.png')             
andesine.reflectance.hi = andesine.reflectance.hi[andesine.reflectance.hi$x > 3,]
andesine.reflectance = rbind(andesine.reflectance.lo, andesine.reflectance.hi)
write.csv(andesine.reflectance, 'usgs/andesine.reflectance.csv')
wavelengths('andesine', 285.0, format(andesine.reflectance$x, digits=5, width=7), format(andesine.reflectance$y, digits=3, width=7) )

# augite.reflectance = digitize::digitize('usgs/splib07a_Augite_NMNH120049_BECKb_AREF.gif.png')                   
# augite.reflectance = digitize::digitize('usgs/splib07a_Augite_WS592_Pyroxene_BECKb_AREF.gif.png')
# augite.reflectance = digitize::digitize('usgs/splib07a_Augite_WS592_Pyroxene_NIC4cc_RREF.gif.png')
# augite.reflectance = digitize::digitize('usgs/splib07a_Augite_WS592_Pyroxene_NIC4cc_RREF_wavenumber.gif.png')
# augite.reflectance = digitize::digitize('usgs/splib07a_Augite_WS588_Pyroxene_NIC4cbb_RREF_wavenumber.gif.png')
augite.reflectance.lo = digitize::digitize('usgs/splib07a_Augite_WS588_Pyroxene_BECKb_AREF.gif.png')               
augite.reflectance.hi = digitize::digitize('usgs/splib07a_Augite_WS588_Pyroxene_NIC4cbb_RREF.gif.png')             
augite.reflectance = rbind(augite.reflectance.lo, augite.reflectance.hi)
write.csv(augite.reflectance, 'usgs/augite.reflectance.csv')
wavelengths('augite', 400.0, format(augite.reflectance$x, digits=5, width=7), format(augite.reflectance$y, digits=3, width=7) )

forsterite.reflectance = digitize::digitize('usgs/splib07a_Forsterite_REE_AZ-01_ASU_ASDFRa_AREF.gif.png')
write.csv(forsterite.reflectance, 'usgs/forsterite.reflectance.csv')
# # grain size is not known, we must use one of the Olivine samples instead
# wavelengths('forsterite', 400.0, format(forsterite.reflectance$x, digits=5, width=7), format(forsterite.reflectance$y, digits=3, width=7) )

# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.a_Fo89_165um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.a_Fo89_165um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.a_Fo89_165um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.b_Fo89_115um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.b_Fo89_115um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.b_Fo89_115um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.c_Fo89_70um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.c_Fo89_70um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.c_Fo89_70um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.d_Fo89_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.d_Fo89_lt60um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.d_Fo89_lt60um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.e_Fo89_lt30um_ASDNGb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.e_Fo89_lt30um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS70.e_Fo89_lt30um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.a_Fo91_65um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.a_Fo91_65um_NIC4cb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.a_Fo91_65um_NIC4cb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.b_Fo91_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.b_Fo91_lt60um_NIC4cb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.b_Fo91_lt60um_NIC4cb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.b_Fo91_lt60um_NIC4cbu_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_GDS71.b_Fo91_lt60um_NIC4cbu_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS285.2B_Fo80_ASDFRb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS285.4B_Fo80_ASDFRb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS285.4B_Fo80_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS285.4B_Fo80_NIC4bb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS285.4B_Fo80_NIC4bb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS420.1B_ASDFRc_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS420.2B_ASDFRc_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS420.3B_ASDFRc_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS420.3B_NIC4cbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_HS420.4B_ASDFRc_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3005_Fo11_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3005_Fo11_lt60um_NIC4bb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3005_Fo11_lt60um_NIC4bb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3054_Fo66_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3054_Fo66_lt60um_NIC4cbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3054_Fo66_lt60um_NIC4cbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3188_Fo51_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3188_Fo51_lt60um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3188_Fo51_lt60um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3189_Fo60_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3189_Fo60_lt60um_NIC4bbu_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3189_Fo60_lt60um_NIC4bbu_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3291_Fo29_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3291_Fo29_lt60um_NIC4bb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3291_Fo29_lt60um_NIC4bb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3377_Fo18_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3377_Fo18_lt60um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI3377_Fo18_lt60um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI4143_Fo41_lt60um_BECKb_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI4143_Fo41_lt60um_NIC4bbb_RREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_KI4143_Fo41_lt60um_NIC4bbb_RREF_wavenumber.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_NMNH137044.a_160u_BECKa_AREF.gif.png')
# olivine.reflectance = digitize::digitize('usgs/splib07a_Olivine_NMNH137044.b_lt74um_BECKa_AREF.gif.png')
olivine.reflectance.lo = digitize::digitize('usgs/splib07a_Olivine_HS420.3B_BECKc_AREF.gif.png')
olivine.reflectance.hi = digitize::digitize('usgs/splib07a_Olivine_HS420.3B_NIC4cbb_RREF.gif.png')       
olivine.reflectance.hi = olivine.reflectance.hi[olivine.reflectance.hi$x > 3,]
olivine.reflectance = rbind(olivine.reflectance.lo, olivine.reflectance.hi)
write.csv(olivine.reflectance, 'usgs/olivine.reflectance.csv')
wavelengths('olivine', mean(c(74,250)), format(olivine.reflectance$x, digits=5, width=7), format(olivine.reflectance$y, digits=3, width=7) )


# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_GDS134_ASDFRb_AREF.gif.png')                   
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_WS219_(limonite)_BECKc_AREF.gif.png')    
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_GDS134_NIC4bbu_RREF.gif.png')                  
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_WS220_BECKc_AREF.gif.png')    
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_GDS134_NIC4bbu_RREF_wavenumber.gif.png')       
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_WS220_NIC4cd_RREF.gif.png')    
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_HS36.3_BECKb_AREF.gif.png')                    
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_WS220_NIC4cd_RREF_wavenumber.gif.png')    
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_HS36.3B_NIC4bcb_RREF.gif.png')                 
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_WS222_Coarse_Gr._BECKa_AREF.gif.png')    
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_HS36.3B_NIC4bcb_RREF_wavenumber.gif.png')      
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_MPCMA2-B_FineGr_adj_BECKb_AREF.gif.png')       
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_MPCMA2-C_M-Crsgrad2_BECKb_AREF.gif.png')       
# goethite.reflectance = digitize::digitize('usgs/splib07a_Goethite_WS222_Medium_Gr._NIC4aau_RREF_wavenumber.gif.png')    
goethite.reflectance.lo = digitize::digitize('usgs/splib07a_Goethite_WS222_Medium_Gr._BECKa_AREF.gif.png')    
goethite.reflectance.hi = digitize::digitize('usgs/splib07a_Goethite_WS222_Medium_Gr._NIC4aau_RREF.gif.png')    
goethite.reflectance.hi = goethite.reflectance.hi[goethite.reflectance.hi$x > 3,]
goethite.reflectance = rbind(goethite.reflectance.lo, goethite.reflectance.hi)
write.csv(goethite.reflectance, 'usgs/goethite.reflectance.csv')
wavelengths('goethite', 125, format(goethite.reflectance$x, digits=5, width=7), format(goethite.reflectance$y, digits=3, width=7) )


# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_GDS483.c_30-60um_ASDFRc_AREF.gif.png')              
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_HS35.3_NIC4cc_RREF_wavenumber.gif.png')
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_GDS483.c_30-60um_NIC4cdb_RREF.gif.png')             
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_S142-1_BECKc_AREF.gif.png')
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_GDS483.c_30-60um_NIC4cdb_RREF_wavenumber.gif.png')  
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_S26-8_BECKc_AREF.gif.png')
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_S29-4_BECKc_AREF.gif.png')
# pyrite.reflectance = digitize::digitize('usgs/splib07a_Pyrite_S30_BECKc_AREF.gif.png')
pyrite.reflectance.lo = digitize::digitize('usgs/splib07a_Pyrite_HS35.3_BECKb_AREF.gif.png')                         
pyrite.reflectance.hi = digitize::digitize('usgs/splib07a_Pyrite_HS35.3_NIC4cc_RREF.gif.png')                        
pyrite.reflectance.hi = pyrite.reflectance.hi[pyrite.reflectance.hi$x > 3,]
pyrite.reflectance = rbind(pyrite.reflectance.lo, pyrite.reflectance.hi)
write.csv(pyrite.reflectance, 'usgs/pyrite.reflectance.csv')
wavelengths('pyrite', mean(c(74,250)), format(pyrite.reflectance$x, digits=5, width=7), format(pyrite.reflectance$y, digits=3, width=7) )


# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.e_20-30um_NIC4bbb_RREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.e_20-30um_NIC4bbb_RREF_wavenumber.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.a_150-250u_NIC4bbb_RREF_wavenumber.gif.png')   
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.f_10-20um_BECKb_AREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.b_104-150u_BECKb_AREF.gif.png')                
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.f_10-20um_NIC4bbb_RREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.b_104-150u_NIC4bbb_RREF.gif.png')              
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.f_10-20um_NIC4bbb_RREF_wavenumber.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.b_104-150u_NIC4bbb_RREF_wavenumber.gif.png')   
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.g_lt10um_BECKb_AREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.c_60-104um_BECKb_AREF.gif.png')                
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.g_lt10um_NIC4dcc_RREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.c_60-104um_NIC4bb_RREF.gif.png')               
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.g_lt10um_NIC4dcc_RREF_wavenumber.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.c_60-104um_NIC4bb_RREF_wavenumber.gif.png')    
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_HS45.3_ASDFRb_AREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.d_30-45um_BECKb_AREF.gif.png')                 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_HS45.3_BECKb_AREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.d_30-45um_NIC4bbb_RREF.gif.png')               
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_HS45.3_NIC4bcu_RREF.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.d_30-45um_NIC4bbb_RREF_wavenumber.gif.png')    
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_HS45.3_NIC4bcu_RREF_wavenumber.gif.png') 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_GDS69.e_20-30um_BECKb_AREF.gif.png')                 
# hematite.reflectance = digitize::digitize('usgs/splib07a_Hematite_WS161_BECKb_AREF.gif.png') 
hematite.reflectance.lo = digitize::digitize('usgs/splib07a_Hematite_GDS69.a_150-250u_BECKb_AREF.gif.png')
hematite.reflectance.hi = digitize::digitize('usgs/splib07a_Hematite_GDS69.a_150-250u_NIC4bbb_RREF.gif.png')
hematite.reflectance.hi = hematite.reflectance.hi[hematite.reflectance.hi$x > 3,]
hematite.reflectance = rbind(hematite.reflectance.lo, hematite.reflectance.hi)
write.csv(hematite.reflectance, 'usgs/hematite.reflectance.csv')
wavelengths('hematite', 233, format(hematite.reflectance$x, digits=5, width=7), format(hematite.reflectance$y, digits=3, width=7) )


# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS195.1B_ASDFRb_AREF.gif.png')              
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.1B_ASDFRb_AREF.gif.png')
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS195.2B_ASDFRb_AREF.gif.png')              
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.2B_ASDFRb_AREF.gif.png')
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS195.3B_ASDFRb_AREF.gif.png')              
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.3B_ASDFRb_AREF.gif.png')
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.3B_BECKb_AREF.gif.png')
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.3B_NIC4bbu_RREF.gif.png')
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS195.3B_NIC4bbu_RREF_wavenumber.gif.png')  
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.3B_NIC4bbu_RREF_wavenumber.gif.png')
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS195.4B_ASDFRb_AREF.gif.png')              
# magnetite.reflectance = digitize::digitize('usgs/splib07a_Magnetite_HS78.4B_ASDFRb_AREF.gif.png')
magnetite.reflectance.lo = digitize::digitize('usgs/splib07a_Magnetite_HS195.3B_BECKb_AREF.gif.png')               
magnetite.reflectance.hi = digitize::digitize('usgs/splib07a_Magnetite_HS195.3B_NIC4bbu_RREF.gif.png')             
magnetite.reflectance.hi = magnetite.reflectance.hi[magnetite.reflectance.hi$x > 3,]
magnetite.reflectance = rbind(magnetite.reflectance.lo, magnetite.reflectance.hi)
write.csv(magnetite.reflectance, 'usgs/magnetite.reflectance.csv')
wavelengths('magnetite', 233, format(magnetite.reflectance$x, digits=5, width=7), format(magnetite.reflectance$y, digits=3, width=7) )


# chalcopyrite.reflectance = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.1B_ASDFRb_AREF.gif.png')   
# chalcopyrite.reflectance = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.2B_ASDFRb_AREF.gif.png')   
# chalcopyrite.reflectance = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.3B_NIC4bc_RREF_wavenumber.gif.png') 
# chalcopyrite.reflectance = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.3B_ASDFRb_AREF.gif.png')   
# chalcopyrite.reflectance = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.4B_ASDFRb_AREF.gif.png') 
# chalcopyrite.reflectance = digitize::digitize('usgs/splib07a_Chalcopyrite_S26-36_BECKb_AREF.gif.png') 
chalcopyrite.reflectance.lo = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.3B_BECKb_AREF.gif.png')    
chalcopyrite.reflectance.hi = digitize::digitize('usgs/splib07a_Chalcopyrite_HS431.3B_NIC4bc_RREF.gif.png') 
chalcopyrite.reflectance.hi = chalcopyrite.reflectance.hi[chalcopyrite.reflectance.hi$x > 3,]
chalcopyrite.reflectance = rbind(chalcopyrite.reflectance.lo, chalcopyrite.reflectance.hi)
write.csv(chalcopyrite.reflectance, 'usgs/chalcopyrite.reflectance.csv')
wavelengths('chalcopyrite', mean(c(74,250)), format(chalcopyrite.reflectance$x, digits=5, width=7), format(chalcopyrite.reflectance$y, digits=3, width=7) )




# anorthite.reflectance = digitize::digitize('usgs/splib07a_Anorthite_GDS28_Syn_lt74um_NIC4aa_RREF_wavenumber.gif.png')  
# anorthite.reflectance = digitize::digitize('usgs/splib07a_Anorthite_HS349.3B_Plagio_NIC4ccc_RREF_wavenumber.gif.png')  
anorthite.reflectance = digitize::digitize('usgs/splib07a_Anorthite_HS201.3B_Plagio_NIC4aa_RREF_wavenumber.gif.png')   

# pyrimidine.reflectance = digitize::digitize('usgs/splib07a_Pyrimidine_SA-131695_82K_NIC4aa_RREF_wavenumber.gif.png')
pyrimidine.reflectance = digitize::digitize('usgs/splib07a_Pyrimidine_SA-131695_NIC4aa_RREF_wavenumber.gif.png')

fluorapatite.reflectance = digitize::digitize('usgs/splib07a_Fluorapatite_WS416_NIC4aa_RREF_wavenumber.gif.png')     

# print.vectors('ethane', format(ethane.reflectance$x, digits=5, width=7), format(ethane.reflectance$y, digits=3, width=7) )
# print.vectors('pyrimidine', format(pyrimidine.reflectance$x, digits=5, width=7), format(pyrimidine.reflectance$y, digits=3, width=7) )

# absorption coefficient (alpha) is related to extinction coefficient (k) and wavenumber (x) by alpha=4πnkx


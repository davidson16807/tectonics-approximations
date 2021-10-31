# WARNING: use the subplot to digitize this line!
lo065=digitize::digitize('loPr.png', x1=0,x2=0.1,y1=0.9,y2=1)
lo070=digitize::digitize('loPr.png', x1=0,x2=0.1,y1=0.9,y2=1)
lo075=digitize::digitize('loPr.png', x1=0,x2=0.1,y1=0.9,y2=1)

# Now do the rest of the plot the usual way
lo060=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo080=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo085=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo090=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo095=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo100=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo105=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo110=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo115=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
# lo120=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo130=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
# lo140=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo160=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)
lo300=digitize::digitize('loPr.png', x1=0,x2=1,y1=0.4,y2=1)


mid100=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid105=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid110=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid115=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid120=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid130=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid140=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
# mid150=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid160=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
# mid180=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid200=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid250=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
# mid350=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)
mid500=digitize::digitize('midPr.png', x1=0,x2=7,y1=0.3,y2=1.1)


hi100=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi120=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi500=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)

hi140=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi160=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi200=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi250=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi300=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi1000=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)
hi1500=digitize::digitize('hiPr.png', x1=0,x2=40,y1=1,y2=4)


lo060$T = 0.60
lo065$T = 0.65
lo070$T = 0.70
lo075$T = 0.75
lo080$T = 0.80
lo085$T = 0.85
lo090$T = 0.90
lo095$T = 0.95
lo100$T = 1.00
lo105$T = 1.05
lo110$T = 1.10
lo115$T = 1.15
# lo120$T = 1.20
lo130$T = 1.30
# lo140$T = 1.40
lo160$T = 1.60
lo300$T = 3.00
mid100$T = 1.00
mid105$T = 1.05
mid110$T = 1.10
mid115$T = 1.15
mid120$T = 1.20
mid130$T = 1.30
mid140$T = 1.40
# mid150$T = 1.50
mid160$T = 1.60
# mid180$T = 1.80
mid200$T = 2.00
mid250$T = 2.50
# mid350$T = 3.50
mid500$T = 5.00
hi100$T = 1.00
hi120$T = 1.20
hi500$T = 5.00
hi140$T = 1.40
hi160$T = 1.60
hi200$T = 2.00
hi250$T = 2.50
hi300$T = 3.00
hi1000$T = 10.00
hi1500$T = 15.00


pTZ=rbind(
	lo060, 
	lo065, 
	lo070, 
	lo075, 
	lo080, 
	lo085, 
	lo090, 
	lo095, 
	lo100, 
	lo105, 
	# lo110, 
	lo115, 
	# lo120, 
	lo130, 
	# lo140, 
	lo160, 
	lo300, 
	mid100, 
	mid105, 
	mid110, 
	mid115, 
	mid120, 
	mid130, 
	mid140, 
	# mid150, 
	mid160, 
	# mid180, 
	mid200, 
	mid250, 
	# mid350, 
	mid500, 
	hi100, 
	hi120, 
	hi500, 
	hi140, 
	hi160, 
	hi200, 
	hi250, 
	hi300, 
	hi1000, 
	hi1500 
)
names(pTZ) = c('p','Z','T')

write.csv(pTZ, file='pTZ.csv')

# Save multiple data objects to working directory
save.image("digitized.RData")
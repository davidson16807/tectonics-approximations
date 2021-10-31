# This script provides functions that tries linear regression on a data set multiple times,
# each time expressing the variables within the data set under different scales (e.g. linear, log, inverse)
# It runs through all possibilities from a given list of scales and returns the combinations 
# that result in the best R^2 value

cartesian.power.list = function(x, n) {
	df = do.call("expand.grid", rep(list(seq(x)), n))
	result = lapply(seq(nrow(df)), function(i) { 
		lapply(seq(ncol(df)), function(j) {
			x[[df[i,j]]]
		}) 
	})
	return(result)
}

test.scales = function(model.list, scales, best.count=3) {
	tests = cartesian.power.list(scales, length(model.list))
	r2s = sapply(tests, function(test) {
		df = do.call(data.frame, lapply(seq(model.list), function(i) {  test[[i]](model.list[[i]])  }) )
		return(summary(lm(df))$r.squared)
	})
	best = tests[order(-r2s)]
	r2s = r2s[order(-r2s)]
	print(r2s)
	for (i in seq(best.count)) {
		print(r2s[i])
		print(best[i])
	}
}

model.list = list(Z, p/T, T)
test.scales(model.list, scales)

scales = list(
	function(x) {x**2},
	function(x) {x**1},
	function(x) {x**0.5},
	function(x) {x**-0.5},
	function(x) {x**-1},
	function(x) {x**-2},

	function(x) {log(x)**2},
	function(x) {log(x)**1},
	function(x) {-log(x)**1},
	function(x) {-log(x)**2}
)
mask = 0.01<I&I<0.99; test.scales(list(I[mask], p[mask], T[mask]), scales)
mask = 0.01<I&I<0.99; plot_ly(x=log(p)[mask], y=(T**-0.5)[mask], z=(I**0.5)[mask])

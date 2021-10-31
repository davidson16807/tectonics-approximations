
You might remember something called the "ideal gas law" from taking high school chemistry class. 

PV = RT

The reason it's called "ideal" is that doesn't always work in reality. Most high school teachers leave it at that but if you take a course in college you might find there's more to the story. Chemists actually know the amount by which a typical gas deviates from the ideal gas law and can predict it to some extent. This is done using something called the the "compressibility factor", which expresses just how far off the mark we are when using the ideal gas law:

Z = PV/RT

In other words, the compressibility factor is the ratio between what we predict ("RT") and what we get ("PV")

The neat thing about the compression factor is this: some guy back in the 1800s discovered it's the same for all gases, just as long as you express the pressure and temperature relative to the pressure and temperature for the critical point of the gas we're dealing with.

Provided we're working in the right units, most gases can be described using values that are measured empirically and compiled for us into a single convenient chart. The most popular of these charts was compiled back in the 1950s by two gentlemen, Nelson and Obert. (This chart is still not quite accurate, most of the remaining error can be corrected by considering what's known as the "assymetric factor" of the molecule, but already we're exceeding the precision needed for the application)

A chart is nice for most researchers but it isn't going to work for me. Sure, I can happily sample data points from that chart and interpolate between them but interpolation is slow. I have a lot of pressures and temperatures to work with and I need results back in less then 1/60 of a second.

So I need some function to approximate it. No problem. There are plenty of approximations for the compressibility factor, some of which are even simple enough to run in my model. However the simplest approximations all have the nasty tendency of limiting their applicability to a particular region in the chart. Remember: I'm trying to simulate whatever I can throw at it. 

So I thought I'd try my hand at writing my own approximation. "No problem", I said. We're talking about a function with only two parameters and one output. I've already done some work in the past to approximate functions for that exact same scenario. How hard could it be?

So I digitized the charts using a library in R that I have some familiarity with, then I adapted the work I did in the past for estimating multivariate Pade approximates.

Pade approximates are beautiful things and I think they're criminally underused. Most people are more familiar with the equations from which they are derived, the polynomial. 

A polynomial can be thought of as a type of formula which captures the full expressivity of addition, subtraction, and multiplication. Any formula that is composed strictly of additions, subtractions, and multiplications can be expressed using a polynomial. If any of you are familiar with abstract algebra, a polynomial can be used to express members of a ring. 

Pade approximates can be thought of as extensions to the polynomial: in addition to addition, subtraction, and multiplication, they also handle division. Despite the added support for division they are nevertheless simple. A pade approximate is nothing more than one polynomial divided by another. Go ahead and try it: write some formula from basic arithmetic and see how it reduces to a ratio of polynomials. If you're the kind who like writing proofs for things, I happened to find the proof for this behavior oddly satisfying. Or, if you like abstract algebra, you can prove to yourself how a pade approximate can be used to express members of a field.

Pade approximates are typically used as an improvement over the Taylor series, since there is a algorithm that converts between them. Surprisingly, the algorithm resembles the algorithm used by Euclid to find ... However, I don't use this algorithm as often as I probably should, since I find the implementation difficult and haven't found a problem that's worth implementing it for. 

More recently however, Pade approximates have found use in to approximate functions even when the Taylor series is not known. This approach was pioneered by Schlick (1994) while writing approximations for functions that are useful to physically based rendering. Given certain aspects of a function, such as the behavior as it approaches infinity, the behavior near 0, and whether or not there are any poles or asymptotes, it's possibly to predict what the pade approximate looks like without even finding the Taylor series or even knowing the derivative. For instance: 

* If a function has an asymptote for a certain value of x, we know the denominator for the pade approximate is zero there
* If a function doesn't have an asymptote anywhere, we know the denominator for the pade approximate must always be either all positive or all negative. The latter case is equivalent to a positive denominator with a negated numerator.
* If a function is symmetric, we know the Pade approximate only has nonzero coefficients for terms with an even number for a degree
* If we know a function behaves linearly as it approaches infinity, we know the largest term in the numerator is one degree higher than the largest term in the divisor. Similar statements can be said for other relationships using L'Hopital's rule.
* If we know a function returns a large output for small input, then the lower degree coefficients of the numerator must be large or the lower degree coefficients of the denominator must be small. Similar statements can be said for large output for large input, small output for large input, or small output for small input. 

Once you know enough about how your function behaves, it becomes easy to weed out the unecessary terms of a pade approximate. From there, you can simply use your favorite optimization algorithm to find the coefficients that give the best approximation of your function. 

But wait! There's more. What happens if you want to approximate a function with more than one parameter? Simple! We already mentioned how any formula can be written as a pade approximate, provided it is only composed of the basic arithmetic operators (➕➖✖➗) A univariate pade approximate starts out by applying these number to combinations of real numbers in addition to a single variable, such as x. What happens if we have more than one variable?

We simply start out with a set of real numbers, plus our first variable, plus our second variable.

I've figured out it's rather straight forward to approximate 2 variable functions with a pade approximate of degree 4 using genetic algorithms. Gradient descent may work for some approximations of lower degree. A degree 4 pade approximate has 30 coefficients, so any approximate with a comparable number of coefficients could also be done.

When I apply my approach to the compressibility factor, I find I can get a mean absolute error of no more than 6%. This puts it in the same ranks as the Papay approximation, mentioned [here](). 
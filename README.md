# Two Series Factorial

This repository will contain information with respect to a function that I created that deals with generating designs for factorial designs dealing with factors, each with two levels. To recall the characteristics of such a design, some brief notes follow below:
*Note that the following content has been written in markdown and has not been rendered to correctly display the intended appearance of the notebook.

# Two-series Designs

Use $2^k$ treamtments; k factors with 2 levels each. We can generate confounding factors if we do not have enough degrees of freedom for estimating of pure error.

The choice of confounding is demonstrated as a $2^{k-q}$ design, where we would like to confound q factors. We confound these variables with blocks; we are not interested in testing blocks. The number of blocks generated by the confounders is $2^q$. 

# Fractional Factorials

In the same vein, we select a subset of a $2^k$ design, by choosing a k - p fractional design, which is denoted as $2^{k-p} = 2x2..x2 k$ p times (which is why a $2^{k-p} fractional design generates a 1/2p fractional design).

In the same way, we need p defining contrasts to choose to be confounded with blocks. The number of blocks generated by this choice of fractional design is the same as if we were choosing q in the two series design; there are $2^p$ blocks of size $2^{k-p}$.

In this case, we select generators $W_1,...,W_p$, which are all aliased with a generic variable named I. 

# Base Factorials and Design Construction.

If we are given the set I, then we can construct base factorials using a sequence of factors that do not all appear concurrenttly in any single element in I. If the concatenation of the factors forms an element in the set I, then we have done it incorrectly.

Alternatively, if we are not given the set I, then we must generate the set I for ourselves. In the usual construction, define k main effects and enumerate all the possible pairs of combinations (standard order) up to $2^{k-p}$ terms. After that, we want to define contrasts. If we are not given a resolution, that we have some maneuverability in our choice of elements of I; although we should try to select elements of I to have high order interaction terms. If we are given a resolution, say resolution(I) = t, then we should define the remaining p main effects as a t-1 order interaction (more on this later).

# Alias of effects of I

To find all possible aliases of the generator I, we need to look at all the elements of the sets belonging to I, and take the combinatorically using addition modulo 2 over the exponents, where we take combinations of size 1, 2, ... , resolution(I) at a time. The defintinion of resolution(I) is the minimal length an element in the set I.

# Alias of Any Effect not I

To find all the possible aliases of a main effect (of which there are k-p many of them), we take the main effect and use addition modulo 2 with all the elements of I. For instance, for main effect of A, we take A*I to find all effects aliased with A.

To define a design with a specified resolution, we need to devise elements of I such that they have length 4. In a table construction that begins with enumerating all the high-low combinations of k-p factors, we devise a way to express the remaining p factors omitted from the main effects by devising a scheme that confounds them that obeys the specified resolution of the design. For instance, if we want a resolution IV (4) design, and we have 7 factors, choosing a subset of 7 - 2 = 5 factors, then we look at main effects A, B, C, D , and E; we need to devise a way to generate the values for F and G. For a resolution 4 design, we can choose F = ABC (because if we take the inverse of F, we can get I = ABCF, which is certainly resolution 4); similarly, we can take G = DEB or some other variation.

An important corollary is if we $begin$ with an aliased set I, we should not consider main effects that utilize every single letter of an element of I; e.g. if I = ABC,...,DEF, we should not generate contrasts that use main effects A,B,C, or D,E,F concurrently.

# Blocks in Fractional Factorials.

We have seen in two-series designs that if we used a confounder, it would generate blocks which would be confounded with the confounder. In partiuclar, if we used q confounders, we would generate $2^q$ blocks, whose elements are treatment terms (high-low factor combinations) deposited into a block based on the parity of the treatment term relative to the defining contrast/confoiunder.

Similarly, given p terms that are aliased, we will generate $2^p$ blocks of size $2^{k-p}$. The corollary is that the treatment terms are aliased, so that there is some mapping of treatment terms to unaliased treatment terms in the model. The specification of which terms belong in which block depends on the parity of the elements of I; they are similarly denoted as -1/+1 (although note that we do not know which is even or odd, it is enough to generate the pairs -/-, -/+, +/-, and +/+) to determine the blocks. So for generator I = -$W_1$, we would deposit half the terms into a block where we have +$W_1$, and the other half the the terms into a block where we have -$W_1$. This generalizes to the use of more defining contrasts (say +DFG/+CEG, +DFG/-CEG, -DFG/+CEG, -CFG/-CEG) for a quarter fractional design.

# Determination of large effects

We have to note that the cost of using fractional designs is a penalty towards our ability to determine which effects are truly significant; larger fractional designs confer more opportunities for variables to become aliased with each other. That means if we have large effects in A, B, C, ..., AB, ...; we have to look at A*I, B*I, .... AB*I to see what other effects have the same name as A, B, C, ... AB, ... 

# De-Aliasing Effects

To determine which effects are TRULY large, we have to devise a new design based on our original design by manipulating the signs of the elements of the set I, More in the text.

# Resolution
As seen before, resolution is defined as the minimal length of an element in the set I. In general, a higher resolution is better than a design with lower resolution because high resolution is associated with the property of confounding with higher order interactions.

# Abberation
An abberation is a count of the number of $k-p+1 ,,, k+p$ way interactions that occur given an aliased set. Generally, designs with lower abberations are better.

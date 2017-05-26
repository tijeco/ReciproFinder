# ReciproFinder
An adorably naiive attempt to find one to one orthologs of N species using all against all blast search

***

## Background

The problem at hand is apparently fairly complicated. The idea is this, given N number of genomes find all one to one ortholog groups. **Why would we want to do that?:** Well it is helpful when making phylogenetic trees, especially if attempting to expand to phylogenomics. This is definitely needed if coalescent aporoaches are needed when you want to deliniate gene trees from the true species tree.
**What are orthologs?:** These are genes shared between common species with a common ancestor, occured due to lineage splitting event. **Why is this hard?** Two reasons: Biologically there are lots of paralogs as well (gene duplications not due to lineage splitting event) which can appear to be orthologous, even though they genuinely aren't and actually have entirely separate evolutionary trajectories. Computationaly, this requires buliding large homology networks and finding groups of reciprocally best homology.

## Approach

**Step 1: all against all blast**

There's unfortunately no whitewashing this step. We have to compare all sequences to each other. That's crazy complex, unavoidable, and only the first step. We literally take our genomes from each species (just coding sequences) as fastas and blast them all agains eachother. To reduce the size of the results file we can limit the E-value to E-5.

**Step2: find reciprocal pairs**

We need to go through the blast file and use bit scores to find orthologous pairs. This should be linear at best, so not too costly (hopefully).

Basic sketch

make hash table
* for i in file

	has gene been put in hash table yet? Put it in there if not along with its partner and accompanying bit score! Other wise check to see if its partner has a better bit score, replacing partners in hash table with best bit score.
	
Now we have who everyone thinks is their bestfriend, now we need to check that it is reciprocated! Go through each added gene, query its partner to see if it matches with itself. Keep track as you go so to only do this once.

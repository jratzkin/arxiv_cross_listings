# arxiv_cross_listings
plot of cross-listings on arxiv, January 2017

This python code accesses metadata on arxiv.org from January 2017, using an API script by Lucas Schwab, available here: 
https://github.com/lukasschwab/arxiv.py . You will need his code (under the file name arxiv.py) to run the python 
code I've written. 

Here is some context for those not familiar with the arxiv. The arxiv, http://arxiv.org, is a preprint server for physics, mathematics, computer science, statistics, quantitative biology, and quantitative finance. Each of those six main categories has within it a number of subcategories (e.g. math.AP, which is the analysis and PDE subcategory in mathematics). When posting a preprint to the arxiv, one must select a primary category, and then one has the option of selecting any number of secondary categories. It is these secondary categories I'm referring to as cross-listing. For instance, one can post a preprint to math.AP and cross-list it on math.DG (differential geometry), to math.DG and math.SP (spectral theory), or not cross-list it at all. 

The output text lists all occuring cross-lists in all subcategories in the month of January 2017. This is a lot of text to absorb, so I made some visuals to go with it. The graphs don't really tell us much, but the bar charts contain some interesting information. For instance, I have bar charts of cross-listings for several subcategories, each of which plots the number of cross-listings in this month. The bar chart of the main categories shows that the ratio of papers with a cross-listing ranges from about .38 to .55, and on average a paper will have about 1.5 cross-listings, among other things. 

Possible improvements: I'm still exploring better ways to visualize the data. 

NB: code takes a while to run because I implemented a 5 second wait time before each request to arxiv.org.

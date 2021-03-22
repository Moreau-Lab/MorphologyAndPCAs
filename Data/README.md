# Explanation of TPS file formatting:
I had a really hard time finding information online about the TPS format that is used for geometric morphometrics, so here is some info for figuring that out:

TPS files consist of blocks of data for each specimen you measure. The measurements section of the block has three columns: the measurement number column, the x coordinate column, and the y coordinate column. These three columns are columns A, B and C for rows 3-16 and 20-33 in the screenshot below. 

Each block of measurements starts with two header cells, one that tells you the number of measurements in each block, formatted as LM=[number of measurements], and one that tells you the scale of the block, formatted as SCALE=[scale value]. These are in cells B1, B2, B18, and B19. For the analysis to work properly, you must have the same number of measurements in every measurement block. 

Finally, the blocks conclude with a cell that gives the specimen ID in the block of measurements above that cell, formatted as ID=[speciment ID]. These are in cells B19 and B34. 

![TPS example](https://github.com/Moreau-Lab/MorphologyAndPCAs/blob/main/Images/TPSexample.png?raw=true)
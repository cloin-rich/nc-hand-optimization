# import libraries
import pandas as pd
import os

# input the path to the folder where you downloaded and extracted the files from the CFIM HAND repository
# CFIM HAND: https://cfim.ornl.gov/data/
# to make this as simple as possible, make sure that the files from CFIM are the only files in the 
# repository, so it should be a collection of folders with different huc6 codes as the folder names

# example path for my desktop
cwd = '/Users/colin/Desktop'

folds = next(os.walk(os.path.join(cwd,'.')))[1]

src_paths = []

# create a list with all of the paths to the src's based on the standard file name
for fold in folds:
    # you only need this if statement if using your desktop
    if fold == '$RECYCLE.BIN':
        continue
    src_path = os.path.join(cwd,fold,'hydrogeo-fulltable-' + fold + '.csv')
    src_paths.append(src_path)

# for each src, read in the file, determine the overlapping COMIDS in two src's,
# drop the duplicate COMIDs from one, and concatenate the two unique src's
for pathi in range(len(src_paths)):
    print('combining src',pathi + 1)
    if pathi == 0: 
        src1 = pd.read_csv(src_paths[pathi])
        src2 = pd.read_csv(src_paths[pathi+1])
    if pathi == 1:
        continue
    if pathi >= 2:
        src1 = src
        src2 = pd.read_csv(src_paths[pathi])   
    catch1 = src1['CatchId'].unique()
    catch2 = src2['CatchId'].unique()
    dupes = np.intersect1d(catch1, catch2)
    src2 = src2[src2.CatchId.isin(np.intersect1d(catch1, catch2)) == False]
    src = pd.concat([src1,src2])

print('finished!')

# nc_src.csv is the name of the output file, but you can change that if you need
src.to_csv(cwd + '/nc_src.csv')
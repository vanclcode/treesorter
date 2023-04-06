# Phylogenetic Tree Analyzer
a seminar program by Adéla Vancl, as part of _NPRG030/Programming I_ at Charles University, Faculty of Mathematics and Physics

## Dependencies
Program requires Python3 with built-in modules sys, os and re

## In-Script Configuration
- `DEFAULT_OUTPUT_FILE` sets default csv output file when unspecified
- `TREE_FILE_EXTENSION` sets file extension to filter when loading directory
- `RUN_QUIETLY` turns off all print functions
- `VERBOSE` turns on some extra informational prints

## Running

Command line interface:

    pta.py [-d DIR][-f FILENAME][-c|-r][-o OUTPUT_FILE] -s SET_DEFINITIONS_STRING -c CRITERIA_DEFINITION_STRING

Example: 

    pta.py -d . -f trees.txt -r -s dino=Dinos-*,karen=Karenia\* -c cany=dino,karen cmax=dino2-,karen3-`

where _trees.txt_ contains

    First-Tree-1000.fasta.tre
    Second-Tree-1150.fasta.tre
    Third-Tree-1200.fasta.tre

## Script Arguments
- **`-d`** sets working directory
- **`-f`** input file containing tree files list
- **`-r`** replace output csv file
- **`-a`** attempt to continue on output file

- **`-c`** string of criteria definitions
- **`-s`** string of set definitions

## String Definitions

### SET_DEFINITIONS_STRING
With format _`name2=First_name,name2=Second_name`_, you can set multiple sets of taxons. 

You can use wildcard **`\*`** anywhere in the name string, such as _`Dinos-\*`_, _`\*-Gracilis\*`_ or _`Kareniaceae-\*-1000\*`_. 

Separate multiple sets by comma without space, names can **alphabetical** characters.

### CRITERIA_DEFINITION_STRING
With format _`crit1=setA,setB, crit2=setA2+,setB5-`_ you can specify criteria to test subtrees for.

To require subtree containing (apart from root taxon) only specific subset of taxons, use comma separated listing of the sets.

To specify minimum and maximum count of set elements in a subtree, append an integer with a + or - after the set. This test if the occurences are more than or equal or less than or equal respectively.

## Output Format
Output is a _csv_ file with values in double quotes.

First line is a header containing column for _root\_taxon_ and columns for each criterion.
At the end of the header row, string of all arguments is added for future reference.

Rows contain first element of tested file (name of the root taxon) and then all maximum bootstraps of specific criterion tests.
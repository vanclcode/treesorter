# Phylogenetic Tree Analyzer
a seminar program by Adéla Vancl, as part of _NPRG030/Programming I_ at Charles University, Faculty of Mathematics and Physics

## Program Specification
- The Program uses defined sets of taxons (which are systematically named _leaves_ of a tree) to **find highest value of an edge** (_bootstrap_), which signifies the probability of the tree being correctly divided by that edge into subtrees of only related _leaves_ – taxons – so that one of the subtrees, apart from specified _seed taxon_, only consists of _taxons_ matching patterns defined by _quantified set operations_ on the aformentioned sets of taxons. The highest _bootstrap_ is systematically saved in an output _CSV file_.

- **CLI of the program** expects arguments of the phylogenetic tree files (_newick_, _nexus_), name of the seed taxon, names of sets and their definition by list of taxon names in combination with wildcards, tested criteria for subtrees defined by quantified set operations, and an output file.

- The program takes each phylogenetic tree file, generates abstract tree structure internally with evaluated edges according to booststrap value, names of the taxons in its leaves. Then it works out defined sets of taxons and for each criterion, it **traverses the tree to find the highest edge value**, which divides the original tree into subtree containing the seed taxon and where all other taxons/leaves match the set criteria. The highest bootstrap value that fits those conditions is saved in the output CSV file in a column respective to the criterion. There can be any number of criteria to search the bootstrap for.

- The conditions/criteria are defined by quantified set operations, which – apart from regular operations – allow for setting minimal and maximal occurence of any specific set members.

- The resulting CSV file contains a row for each tested phylogenetic tree file and columns respective to the criteria set by program arguments.

## Specifikace programu (Czech)
- Program dle definovaných _množin taxonů_ (systematicky pojmenovaných _listů_ stromu) v každém stromu **nalezne nejvyšší ohodnocení hrany** (_bootstrap_), které určuje pravděpodobnost správného rozdělení příbuznosti listů ve dvou **podstromech oddělených touto hranou**, kde jeden ze stromů obsahuje kromě zadaného _seed taxonu_ pouze _taxony_ popsané _kvantifikovanými množinovými operacemi_, a nalezenou hodnotu systematicky uloží do výstupního _CSV souboru_.

- **CLI programu** očekává jako argumenty soubor(y) fylogenetických stromů v závorkovém formátu (_newick_, _nexus_), název seed taxonu, jména množin a výčet (s _wildcardy_) názvů taxonů do množiny spadajících, zkoumané podmínky pro podstromy definované kvantifikovanými množinovými operacemi a výstupní soubor.

- Program pro každý soubor fylogenetického stromu vygeneruje abstraktní strukturu stromu s ohodnocenými hranami a názvy taxonů v listech. Zpracuje zadané výčty množin taxonů a pro každou ze zadaných podmínek **hledá hranu s nejvyšší hodnotou**, která odděluje podstrom obsahující seed taxon a kde všechny ostatní listy spňují řešenou podmínku. Nejvyšší hodnotu ukládá do CSV souboru do příslušného sloupce.

- Podmínky jsou zadány kvantifikovanými množinovými operacemi, kde lze **kromě obvyklých operací definovat i minimální nebo maximální četnost výskytu** prvku množiny.

- Ve výsledném CSV souboru každý **řádek odpovídá jednomu ze souborů stromu** a sloupce odpovídají argumenty zadaným podmínkám.

_This specification has been approved by Mgr. Rudolf Rosa, Ph.D. (seminar program supervisor) and Mgr. Anna M. G. Novák Vanclová (bioinformatics consultant and project owner)_

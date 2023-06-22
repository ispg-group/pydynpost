# Dictonary containing tasks as keys and module name as values
moduleNames = {"population"   : "aimspopulation", 
               "complexity"   : "aimscomplexity",
               "coupling"     : "aimscoupling", 
               "selection"    : "aimsselection",
               "history"      : "aimshistories",
               "internals"    : "aimsexpectationvalues",
               "molpop"       : "aimsmolpop",
               "decoherence"  : "aimsdecoherence",
               "momentum"     : "aimsspawninganalysis"} 

# Dictonary containing tasks as keys and class names as values;
# Important because some modules contain more than one class
methodNames = {"population"   : "statePopulations", 
               "complexity"   : "computationalComplexity",
               "coupling"     : "analyzeCouplings", 
               "selection"    : "selectionConsistency",
               "history"      : "historyStatistics",
               "internals"    : "internals",
               "molpop"       : "molpop",
               "decoherence"  : "decoherenceTimes",
               "momentum"     : "analyzeSpawning"}

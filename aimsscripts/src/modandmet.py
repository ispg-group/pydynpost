# Dictonary containing tasks as keys and module name as values
moduleNames = {"population"   : "population", 
               "complexity"   : "complexity",
               "coupling"     : "coupling", 
               "selection"    : "selection",
               "history"      : "histories",
               "internals"    : "expectationvalues",
               "molpop"       : "expectationvalues",
               "decoherence"  : "decoherence",
               "momentum"     : "spawninganalysis"} 

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

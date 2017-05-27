pairDict = {}
with open("Sim_genomes.fasta.blastall") as f:
    for line in f:
        row = line.strip().split()
        bit_score = float(row[11])
        gene1 = row[0]
        gene2 = row[1]
        taxa1 = gene1.split('_')[0]
        taxa2 = row[1].split('_')[0]
        if taxa1 != taxa2:
            if taxa1 not in pairDict:
                pairDict[taxa1] = {}
            if gene1 not in pairDict[taxa1]:
                pairDict[taxa1][gene1] = {}
            if taxa2 not in pairDict[taxa1][gene1]:
                pairDict[taxa1][gene1][taxa2] = (gene2, bit_score)
            elif pairDict[taxa1][gene1][taxa2][1] < bit_score:
                pairDict[taxa1][gene1][taxa2] = (gene2, bit_score)
    bestieDict = {}
    for taxa in pairDict.keys():
        for gene in pairDict[taxa].keys():
            for genePairTaxa in pairDict[taxa][gene]:
                bestie =  pairDict[taxa][gene][genePairTaxa][0]
                bestieOfBestie = pairDict[genePairTaxa][bestie][taxa][0]
                if gene == bestieOfBestie:
                    if min(bestie,gene) not in bestieDict:
                        bestieDict[min(bestie,gene)] = max(bestie,gene)
                    else:
                        print(bestie, "and",gene,"are BFFs")

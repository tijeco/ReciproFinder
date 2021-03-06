pairDict = {}

genesUsedDict ={}
orthologs = {}
number = 1
def checkPairs(original_dict,internalDict,value):
    for j in internalDict.keys():

        if((max(value,j),min(value,j)) not in orthologs[number]):

            if j in original_dict:
                orthologs[number][(max(value,j),min(value,j))] = True

                checkPairs(original_dict,original_dict[j],j)
            # else:
            orthologs[number][(max(value,j),min(value,j))] = True
            genesUsedDict[j] =True
            genesUsedDict[value] =True
            if j in orthologs[number]["genes"]:
                orthologs[number]["genes"][j] += 1
            else:
                orthologs[number]["genes"][j] = 1
            if value in orthologs[number]["genes"]:
                orthologs[number]["genes"][value] += 1
            else:
                orthologs[number]["genes"][value] = 1

            orthologs[number]["taxa"][j.split("_")[0]] =True
            orthologs[number]["taxa"][value.split("_")[0]] =True
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
                try:
                    bestieOfBestie = pairDict[genePairTaxa][bestie][taxa][0]
                    if gene == bestieOfBestie:
                        if min(bestie,gene) not in bestieDict:
                            bestieDict[min(bestie,gene)] = {max(bestie,gene):True}
                        else:
                            bestieDict[min(bestie,gene)][max(bestie,gene)] = True
                            print(bestie, "and",gene,"are BFFs")
                except:
                    None

    for i in bestieDict.keys():
        orthologs[number] = {}
        orthologs[number]["genes"] ={}
        orthologs[number]["taxa"] ={}
        checkPairs(bestieDict,bestieDict[i],i)
        number+=1
    for i in orthologs.keys():
        # print(len(orthologs[i])-2)
        taxa_completion = len(orthologs[i]["taxa"])/float(len(pairDict))
        edge_count = 0
        taxa_counter = {}
        for j in orthologs[i]["genes"]:
            edge_count+= orthologs[i]["genes"][j]
            taxa_counter[j.split("_")[0]] = True
        taxa_uniq = len(taxa_counter)
        edge_completion =(edge_count/len(orthologs[i]["genes"])/3)
        if edge_completion > 0.65 and taxa_completion > 0.75 and taxa_uniq == len(pairDict) and len(pairDict) == len(orthologs[i]["genes"]):
            print(orthologs[i]["genes"])

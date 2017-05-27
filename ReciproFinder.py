pairDict = {}
with open("sampleBlast.txt") as f:
    for line in f:
        row = line.strip().split()
        bit_score = int(row[2])
        taxa1 = row[0].split('_')[0]
        taxa2 = row[1].split('_')[0]

        # print(row[0] == row[1])
        if taxa1 != taxa2:
            if taxa1 not in pairDict:
                pairDict[taxa1] = {}
                pairDict[taxa1][row[0]] =  (row[1],bit_score)
            else:
                try:
                    if pairDict[taxa1][row[0]][1] < bit_score:
                        pairDict[taxa1][row[0]] =  (row[1],bit_score)
                except:
                    pairDict[taxa1][row[0]] =  (row[1],bit_score)
    print(pairDict)
    for i in pairDict.keys():
        for j in pairDict[i].keys():
            bestie =pairDict[i][j][0]
            bestieTaxa = pairDict[i][j][0].split("_")[0]
            bestieOfBestie = pairDict[bestieTaxa][bestie][0]
            if(j == bestieOfBestie):
                print(j," and ", bestie," are BFFs")

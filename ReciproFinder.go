package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

type orthoGroup struct {
	nodePairs map[Pair]bool
	geneNodes map[string]int
	taxaNodes map[string]bool
}

type Pair struct {
	a, b interface{}
}

func keysInOriginal(originalMap map[string]map[string]bool, subMap map[string]bool, usedGenes map[string]bool, originalGene string, outputStruct orthoGroup) {
	// fmt.Println(originalGene, "==============================")
	fmt.Println("currently working with the following:", originalGene, subMap)
	//NOTE originalMap[originalGene] ======>>>>>> subMap, so yeah man that's a thing NOTE NOTE NOTE

	// usedGenes[originalGene] = true

	// if ok {
	// 	fmt.Println(false)
	// usedGenes[newGene] = true
	for newGene := range subMap {
		fmt.Println("newGene:", newGene, "originalGene:", originalGene)
		bestiePair := []string{newGene, originalGene}
		sort.Strings(bestiePair)

		_, putatitivePairs := outputStruct.nodePairs[Pair{bestiePair[0], bestiePair[1]}]
		if !putatitivePairs {

			// outputStruct.nodePairs[Pair{newGene, originalGene}] = true

			outputStruct.nodePairs[Pair{bestiePair[0], bestiePair[1]}] = true

			taxaList := []string{strings.Split(newGene, "_")[0], strings.Split(originalGene, "_")[0]}
			for thisTaxa := range taxaList {
				outputStruct.taxaNodes[taxaList[thisTaxa]] = true
			} //end taxa loop
			// geneList := []string{newGene, originalGene}
			// fmt.Println(geneList)
			// fmt.Println(outputStruct.geneNodes)
			// for currentGene := range geneList {
			_, ok := outputStruct.geneNodes[newGene]
			if ok {
				outputStruct.geneNodes[newGene]++ //FIXME
				fmt.Println(newGene, outputStruct.geneNodes[newGene])
			} else {
				outputStruct.geneNodes[newGene] = 1
			} //end geneNodes comparison

			_, okay := outputStruct.geneNodes[originalGene]
			if okay {
				outputStruct.geneNodes[originalGene]++ //FIXME
				fmt.Println(originalGene, outputStruct.geneNodes[originalGene])
			} else {
				outputStruct.geneNodes[originalGene] = 1
			} //end geneNodes comparison
			// } //end gene loop
			if putatitiveGene, ok := originalMap[newGene]; ok {
				usedGenes[newGene] = true
				keysInOriginal(originalMap, putatitiveGene, usedGenes, newGene, outputStruct)

			} // end recursion
		} //end putatitivePairs

	} // end submap loop
	fmt.Println("done with function", originalGene)
} // end keysInOriginal

func main() {
	// var numTaxa int
	pairMap := make(map[string]map[string]map[string]Pair)
	// pairMap["A1"] = make(map[string]Pair)
	// pairMap["A1"]["g1"] = Pair{"g2", 33.4}
	if file, err := os.Open("Sim_genomes.fasta.blastall.big"); err == nil {

		// make sure it gets closed
		defer file.Close()

		// create a new scanner and read the file line by line
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			row := strings.Fields(scanner.Text())
			// log.Println(row[0])
			bitScore, _ := strconv.ParseFloat(row[11], 64) //deal with error here, could be important
			gene1 := row[0]
			gene2 := row[1]
			taxa1 := strings.Split(gene1, "_")[0]
			taxa2 := strings.Split(gene2, "_")[0]
			if taxa1 != taxa2 {
				if currentTaxa, ok := pairMap[taxa1]; ok {
					if currentGene, ok := pairMap[taxa1][gene1]; ok {
						// log.Println(currentGene.b.(float64) > bitScore)
						if otherTaxa, ok := pairMap[taxa1][gene1][taxa2]; ok {
							if otherTaxa.b.(float64) < bitScore {
								currentGene[taxa2] = Pair{gene2, bitScore}
							}
						} else {
							currentGene[taxa2] = Pair{gene2, bitScore}
						} // end otherTaxa
					} else {
						currentTaxa[gene1] = make(map[string]Pair)        // initialize submap for this gene
						currentTaxa[gene1][taxa2] = Pair{gene2, bitScore} //the would be bitScore
					} //end currentGene
				} else {
					pairMap[taxa1] = make(map[string]map[string]Pair) // haven't seen this taxa yet, initialize its submap
					pairMap[taxa1][gene1] = make(map[string]Pair)
					pairMap[taxa1][gene1][taxa2] = Pair{gene2, bitScore} // the would be bitScore
				} //end currentTaxa

			} // end taxa comparison
			//NOTE put number of taxa as len(pairDict) here
		} // end file scan loop

	} else {
		log.Fatal(err)
	} //end file open
	bestieMap := make(map[string]map[string]bool)
	for taxa, taxaGene := range pairMap {
		// fmt.Println("taxa:", taxa, taxaGene)
		for gene, geneTaxa := range taxaGene {
			// fmt.Println("gene:", gene)
			for genePairTaxa, genePair := range geneTaxa {
				// fmt.Println("genePairTaxa:", gene, genePairTaxa, genePair)
				bestie := genePair.a.(string)
				// fmt.Println(genePairTaxa, bestie, taxa, gene)
				// fmt.Println(pairMap[genePairTaxa][bestie][taxa].a, gene)

				//NOTE this has potential for key error, consider checking first
				bestieOfBestie, _ := pairMap[genePairTaxa][bestie][taxa].a.(string)

				// fmt.Println(gene == bestieOfBestie)
				if gene == bestieOfBestie {
					biffleFriends := []string{gene, bestie}
					sort.Strings(biffleFriends)
					// fmt.Println("************", biffleFriends)
					if firstGene, ok := bestieMap[biffleFriends[0]]; ok {
						firstGene[biffleFriends[1]] = true

					} else {
						bestieMap[biffleFriends[0]] = make(map[string]bool)
						bestieMap[biffleFriends[0]][biffleFriends[1]] = true
					} //end bestieMap filling

				} // end compare gene to bestieOfBestie

			} //end genePairTaxa loop

		} //end gene loop
	} // end taxa loop
	// fmt.Println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	orthologNumber := 1
	// orthologs := orthoGroup{}
	// orthologs.nodePairs = make(map[Pair]bool)
	// orthologs.geneNodes = make(map[string]int)
	// orthologs.taxaNodes = make(map[string]bool)
	orthologMap := make(map[int]orthoGroup)

	pathTraveled := make(map[string]bool)

	// fmt.Println(bestieMap)
	for gene := range bestieMap {
		// fmt.Println(gene, value, "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

		orthologs := orthoGroup{}
		orthologs.nodePairs = make(map[Pair]bool)
		orthologs.geneNodes = make(map[string]int)
		orthologs.taxaNodes = make(map[string]bool)
		orthologMap[orthologNumber] = orthologs
		// fmt.Println(gene, "***********************************")
		_, ok := pathTraveled[gene]
		if !ok {
			keysInOriginal(bestieMap, bestieMap[gene], pathTraveled, gene, orthologMap[orthologNumber])
			fmt.Println("new row!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		}

		pathTraveled[gene] = true
		orthologNumber++
	}
	fmt.Println("??????????????????????????")
	for gene, value := range bestieMap {
		fmt.Println(gene, value)

	}
	fmt.Println("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", orthologMap)
	for num, value := range orthologMap {
		fmt.Println(num, value.geneNodes)
	}

	strs := []string{"c", "a"}
	sort.Strings(strs)
	fmt.Println("Strings:", strs[0])

} // end main

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
	nodePairs     map[Pair]bool
	taxaGeneNodes map[string]map[string]int
}

type Pair struct {
	a, b interface{}
}

func keysInOriginal(originalMap map[string]map[string]bool, subMap map[string]bool, usedGenes map[string]bool, originalGene string, outputStruct orthoGroup) {
	for newGene := range subMap {
		// fmt.Println("newGene:", newGene, "originalGene:", originalGene)
		bestiePair := []string{newGene, originalGene}
		sort.Strings(bestiePair)

		_, putatitivePairs := outputStruct.nodePairs[Pair{bestiePair[0], bestiePair[1]}]
		if !putatitivePairs {

			outputStruct.nodePairs[Pair{bestiePair[0], bestiePair[1]}] = true

			if putatitiveTaxa, ok := outputStruct.taxaGeneNodes[strings.Split(newGene, "_")[0]]; ok {

				_, putatitiveGene := putatitiveTaxa[newGene]
				if putatitiveGene {
					putatitiveTaxa[newGene]++
				} else {
					putatitiveTaxa[newGene] = 1
				}
			} else {
				outputStruct.taxaGeneNodes[strings.Split(newGene, "_")[0]] = map[string]int{newGene: 1}
			}

			if putatitiveTaxa, ok := outputStruct.taxaGeneNodes[strings.Split(originalGene, "_")[0]]; ok {
				_, putatitiveGene := putatitiveTaxa[originalGene]
				if putatitiveGene {
					putatitiveTaxa[originalGene]++
				} else {
					putatitiveTaxa[originalGene] = 1
				}

			} else {
				// putatitiveTaxa[originalGene] = 1
				outputStruct.taxaGeneNodes[strings.Split(originalGene, "_")[0]] = make(map[string]int) //{originalGene: 1}
				outputStruct.taxaGeneNodes[strings.Split(originalGene, "_")[0]][originalGene] = 1
			}

			if putatitiveGene, ok := originalMap[newGene]; ok {
				usedGenes[newGene] = true
				keysInOriginal(originalMap, putatitiveGene, usedGenes, newGene, outputStruct)

			} // end recursion
		} //end putatitivePairs

	} // end submap loop
} // end keysInOriginal

func main() { 
	var numTaxa int
	pairMap := make(map[string]map[string]map[string]Pair)
	arg := os.Args[1]
	// arg := "Sim_genomes.fasta.blastall"
	if file, err := os.Open(arg); err == nil {
		// NOTE this opens file
		// fmt.Println("opening:", arg)

		// make sure it gets closed
		defer file.Close()

		// create a new scanner and read the file line by line
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			row := strings.Fields(scanner.Text())
			bitScore, _ := strconv.ParseFloat(row[11], 64) //deal with error here, could be important
			gene1 := row[0]
			gene2 := row[1]
			taxa1 := strings.Split(gene1, "_")[0]
			taxa2 := strings.Split(gene2, "_")[0]
			if taxa1 != taxa2 {
				if currentTaxa, ok := pairMap[taxa1]; ok {
					if currentGene, ok := pairMap[taxa1][gene1]; ok {
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
		numTaxa = len(pairMap)

	} else {
		log.Fatal(err)
	} //end file open
	// NOTE this closes file
	// fmt.Println("closing:", arg)
	bestieMap := make(map[string]map[string]bool)
	for taxa, taxaGene := range pairMap {
		for gene, geneTaxa := range taxaGene {
			for genePairTaxa, genePair := range geneTaxa {
				bestie := genePair.a.(string)

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
	orthologNumber := 1

	orthologMap := make(map[int]orthoGroup)

	pathTraveled := make(map[string]bool)

	for gene := range bestieMap {

		orthologs := orthoGroup{}
		orthologs.nodePairs = make(map[Pair]bool)
		orthologs.taxaGeneNodes = make(map[string]map[string]int)
		orthologMap[orthologNumber] = orthologs
		_, ok := pathTraveled[gene]
		if !ok {
			keysInOriginal(bestieMap, bestieMap[gene], pathTraveled, gene, orthologMap[orthologNumber])
		}

		pathTraveled[gene] = true
		orthologNumber++
	}

	for num, value := range orthologMap {

		if len(value.taxaGeneNodes) == numTaxa {
			// NOTE this is alot of stuff to print
			// fmt.Println(num, (float64(len(value.nodePairs))*2.0)/float64(numTaxa) > 0.6*float64(numTaxa-1), value.taxaGeneNodes)
			if (float64(len(value.nodePairs))*2.0)/float64(numTaxa) > 0.6*float64(numTaxa-1) {
				for taxa := range value.taxaGeneNodes {
					genes := value.taxaGeneNodes[taxa]
					var winner string
					var highestNum int
					for gene, connections := range genes {
						if len(genes) == 1 {
							fmt.Println(num, gene)
							// fmt.Println(num, gene, connections)
						} else {
							// fmt.Println("+++++++++++++", genes, gene, taxa)
							if connections > highestNum {
								winner = gene
								highestNum = connections
							}
						}
					}
					if winner != "" {
						fmt.Println(num, winner)
						// fmt.Println(num, winner, highestNum)
					}
				}

			}

		}
	}
} // end main

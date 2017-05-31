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
	node_pairs map[Pair]bool
	gene_nodes map[string]int
	taxa_nodes map[string]bool
}

type Pair struct {
	a, b interface{}
}

func keysInOriginal(originalMap map[string]map[string]bool, subMap map[string]bool, originalGene string, outputStruct orthoGroup) {

	for newGene := range subMap {
		outputStruct.node_pairs[Pair{newGene, originalGene}] = true
		taxaList := []string{strings.Split(newGene, "_")[0], strings.Split(originalGene, "_")[0]}
		for thisTaxa := range taxaList {
			outputStruct.taxa_nodes[taxaList[thisTaxa]] = true
		} //end taxa loop
		geneList := []string{originalGene, newGene}
		for currentGene := range geneList {
			_, ok := outputStruct.gene_nodes[geneList[currentGene]]
			if ok {
				outputStruct.gene_nodes[geneList[currentGene]]++
			} else {
				outputStruct.gene_nodes[geneList[currentGene]] = 1
			} //end gene_nodes comparison
		} //end gene loop
		if putative_gene, ok := originalMap[newGene]; ok {
			keysInOriginal(originalMap, putative_gene, newGene, outputStruct)
		} // end recursion
	} // end submap loop
} // end keysInOriginal

func main() {
	pairMap := make(map[string]map[string]map[string]Pair)
	// pairMap["A1"] = make(map[string]Pair)
	// pairMap["A1"]["g1"] = Pair{"g2", 33.4}
	if file, err := os.Open("Sim_genomes.fasta.blastall"); err == nil {

		// make sure it gets closed
		defer file.Close()

		// create a new scanner and read the file line by line
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			row := strings.Fields(scanner.Text())
			log.Println(row[0])
			bit_score, _ := strconv.ParseFloat(row[11], 64) //deal with error here, could be important
			gene1 := row[0]
			gene2 := row[1]
			taxa1 := strings.Split(gene1, "_")[0]
			taxa2 := strings.Split(gene2, "_")[0]
			if taxa1 != taxa2 {
				if current_taxa, ok := pairMap[taxa1]; ok {
					if current_gene, ok := pairMap[taxa1][gene1]; ok {
						// log.Println(current_gene.b.(float64) > bit_score)
						if other_taxa, ok := pairMap[taxa1][gene1][taxa2]; ok {
							if other_taxa.b.(float64) < bit_score {
								current_gene[taxa2] = Pair{gene2, bit_score}
							}
						} else {
							current_gene[taxa2] = Pair{gene2, bit_score}
						} // end other_taxa
					} else {
						current_taxa[gene1] = make(map[string]Pair)         // initialize submap for this gene
						current_taxa[gene1][taxa2] = Pair{gene2, bit_score} //the would be bit_score
					} //end current_gene
				} else {
					pairMap[taxa1] = make(map[string]map[string]Pair) // haven't seen this taxa yet, initialize its submap
					pairMap[taxa1][gene1] = make(map[string]Pair)
					pairMap[taxa1][gene1][taxa2] = Pair{gene2, bit_score} // the would be bit_score
				} //end current_taxa

				bestieMap := make(map[string]map[string]bool)
				for taxa, taxa_genes := range pairMap {
					fmt.Println("taxa:", taxa, taxa_genes)
					for gene, gene_taxa := range taxa_genes {
						fmt.Println("gene:", gene)
						for genePairTaxa, gene_pair := range gene_taxa {
							fmt.Println("genePairTaxa:", gene, genePairTaxa, gene_pair)
							bestie := gene_pair.a.(string)
							fmt.Println(genePairTaxa, bestie, taxa, gene)
							fmt.Println(pairMap[genePairTaxa][bestie][taxa].a, gene)

							//NOTE this has potential for key error, consider checking first
							bestieOfBestie, _ := pairMap[genePairTaxa][bestie][taxa].a.(string)

							fmt.Println(gene == bestieOfBestie)
							if gene == bestieOfBestie {
								biffleFriends := []string{gene, bestie}
								sort.Strings(biffleFriends)
								fmt.Println("************", biffleFriends)
								if first_gene, ok := bestieMap[biffleFriends[0]]; ok {
									first_gene[biffleFriends[1]] = true

								} else {
									bestieMap[biffleFriends[0]] = make(map[string]bool)
									bestieMap[biffleFriends[0]][biffleFriends[1]] = true
								} //end bestieMap filling

							} // end compare gene to bestieOfBestie

						} //end genePairTaxa loop

					} //end gene loop
				} // end taxa loop
				fmt.Println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
				ortholog_number := 1
				orthologs := orthoGroup{}
				orthologs.node_pairs = make(map[Pair]bool)
				orthologs.gene_nodes = make(map[string]int)
				orthologs.taxa_nodes = make(map[string]bool)
				orthologMap := make(map[int]orthoGroup)
				fmt.Println(bestieMap)
				for gene := range bestieMap {
					orthologMap[ortholog_number] = orthologs
					keysInOriginal(bestieMap, bestieMap[gene], gene, orthologMap[ortholog_number])
					ortholog_number++
				}
				fmt.Println("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", orthologMap)

			} // end taxa comparison
			//NOTE put number of taxa as len(pairDict) here
		} // end file scan loop

	} else {
		log.Fatal(err)
	} //end file open

	strs := []string{"c", "a"}
	sort.Strings(strs)
	fmt.Println("Strings:", strs[0])

} // end main

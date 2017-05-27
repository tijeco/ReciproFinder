
import random
codons = ["CTA","CGT","TGC","CAT"]





def randCodon():
    stops = ["TAA","TAG","TGA"]
    codon = "TAA"
    while codon in stops:
        #if codon in stops:
            codon= random.choice("ATCG")+random.choice("ATCG")+random.choice("ATCG")
            if codon not in stops:
                    return codon

seq=""
for i in range(100):
    seq+=randCodon()
for i in range(4):
    with open("gene"+str(1+i)+".fa","w") as out:
        out.write("ATG"+seq+"TAA\n")

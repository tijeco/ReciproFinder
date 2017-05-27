number = 1
with open("seed_Unaligned.FASTA") as f:
    if line[0] != ">":
        with open("gene"+str(number)+".fa","w") as out:
            out.write(line)
        

    number+=1

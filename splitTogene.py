import os

number = 1
with open("seed_Unaligned.FASTA") as f:
    for line in f:
        if line[0] != ">":
            filename = "gene"+str(number)+"/gene"+str(number)+".fa"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            with open(filename,"w") as out:
                out.write(line)


            number+=1

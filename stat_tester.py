#!/home/yigitaras/anaconda3/bin/python3
import sys , subprocess

ex_name = "RAND_DIST.cpp"
ex_name2 = "GENETIC_DIST.cpp"


n = int(sys.argv[1])
p = int(sys.argv[2])

subprocess.call(['g++' , ex_name,"-o","rand_gen","-O3"])
subprocess.call(["./rand_gen",str(n),str(p)])
subprocess.call(["g++", ex_name2,"-o",'gen_gen',"-O3"])
subprocess.call(["./gen_gen",str(n),str(p)])


f1Name = "random_dist_" + str(n) + "_" + str(p) + ".txt"
f2Name = "genetic_dist_" + str(n) + "_" + str(p) + ".txt"
ct1List = [0]*(1+(n-1)**2)
ct2List = [0]*(1+(n-1)**2)

inFile1 = open(f1Name,"r")
inFile2 = open(f2Name,"r")


with open(f1Name,"r") as f:
    temp = ""
    while True:
        c = f.read(1)
        if not c:
            break
        if c != " ":
            temp += c
        else:
            ct1List[int(temp)] += 1
            temp = ""
            
with open(f2Name,"r") as f:
    temp = ""
    while True:
        c = f.read(1)
        if not c:
            break
        if c != " ":
            temp += c
        else:
            ct2List[int(temp)] += 1
            temp = ""



f = open("stat_test_result.txt", "w")

f.write("Rand\t\t\tGen\n")
for k in range(len(ct1List)):
    f.write("len("+str(k)+")" + ":" + str(ct1List[k]) + "\t\t\t" + "len("+str(k)+")"+ ":" + str(ct2List[k])+"\n")
    f.write("=" * 40 + "\n" )


f.close()

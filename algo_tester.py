#!/home/yigitaras/anaconda3/bin/python3

import sys, subprocess, os

the_exec = "./Greedy_Gen"
the_exec2 = "./Genetic_Gen"
N = sys.argv[1]
P = sys.argv[2]
pop = "50"
cycle = "300"
cut_off = "3"
print("Script Started")
for _ in range(0,300):
    op1 = subprocess.call([the_exec, N , P])
for _ in range(0,300):
    op2 = subprocess.call([the_exec2,N,P,pop,cycle,cut_off])
file_name = "greedy_log_"+str(N)+"_"+str(P)+".txt"
file_name2 = "genetic_log_"+str(N)+"_"+str(P)+".txt"
greedy_beg_avg = 0
greedy_end_av = 0
with open(file_name,"r") as fp:
    greedy_beg_val = 0
    greedy_end_val = 0
    ctr = 0
    for line in fp:
        vals = line.split()
        greedy_beg_val += int( vals[0] )
        greedy_end_val += int( vals[1] )
        ctr += 1
    greedy_beg_avg = greedy_beg_val/ctr
    greedy_end_avg = greedy_end_val/ctr
genetic_beg_avg = 0
genetic_end_avg = 0
with open(file_name2,"r") as fp:
    genetic_beg_val = 0
    genetic_end_val = 0
    ctr = 0
    for line in fp:
        vals = line.split()
        genetic_beg_val += int( vals[0] )
        genetic_end_val += int( vals[1] )
        ctr += 1
    genetic_beg_avg = genetic_beg_val/ctr
    genetic_end_avg = genetic_end_val/ctr
greedy_beg_avg = round(greedy_beg_avg,2)
greedy_end_avg = round(greedy_end_avg,2)
genetic_beg_avg = round(genetic_beg_avg,2)
genetic_end_avg =round(genetic_end_avg,2)

test_res_file = "test_results.txt"

check = os.path.exists(test_res_file)
if check is False:
    f_header = open(test_res_file,"w")
    f_header.write("The Test Cases For Genetic And Greedy Algorithms\n")
    f_header.write("\n")
    f_header.write("The tests for Genetic Algorithm has been done with 25 population, 300 cycles and a 1/3 cut off rate.Average is taken from 300 iterations.\n")
    f_header.write("\n")
    f_header.close()

with open(test_res_file,"a") as fp:
    fp.write("For Nodes:{} and Alphabet Size:{}. ".format(N,P))
    fp.write("The AVG of BEG and END for Greedy and Genetic: {} .{} , {} {}\n".format(greedy_beg_avg,greedy_end_avg,
                                                                                   genetic_beg_avg,genetic_end_avg))


print("Done")

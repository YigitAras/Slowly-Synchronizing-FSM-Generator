#include <iostream>
#include <queue>
#include <cmath>
#include <climits>
#include <string.h>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <random>
#include <stdlib.h>

#define DEBUG

#define Id(s1, s2) ((s1 > s2)?(((s1 * (s1 + 1))/2) + s2):(((s2 * (s2 + 1))/2) + s1)) //this is how we compute the ids
#define IdNaive(s1, s2, N) ((s1 > s2)?((s1 * N) + s2):((s2 * N) + s1)) //this is how we compute the ids
#define s1fromId(id) ((int)(sqrt((2.0 * id) +1.0) - 0.5));
#define s2fromId(id, s1) (id - ((s1 * (s1 + 1))/2));

using namespace std;


#define MUTATE_SELF 0;
#define MUTATE_INDEG 1;
#define THREE_CROSS 2;
#define NODE_WISE 3;
#define COMP_RAND 4;
#define RAND_MUT 5;
#define ALPH_WISE 6;

struct PNode{
    int letter;
    PNode* next;
    PNode(int _letter, PNode* _next) : letter(_letter), next(_next) {}
};

template <typename T>
struct DistID { //ege: to store pairs with their distances in min heap
    T id;
    T dist;
    T label;
    DistID(T i, T d, T l) : id(i), dist(d), label(l) {}

    bool operator<(const DistID& rhs) const {
        if (dist > rhs.dist)
            return true;
        else if (dist < rhs.dist)
            return false;
        else {
            if (label > rhs.label)
                return true;
            else
                return false;
        }//return dist > rhs.dist;
    }
};

bool syncCheck(int* a, int* iap, int* ia, int N, int P) {
    int noOfPair = (N * (N + 1)) / 2;

    int* distance = new int[noOfPair];
    int* next = new int[noOfPair];
    int* letter = new int[noOfPair];
    int* que = new int[noOfPair];

    for (int i = 0; i < noOfPair; i++) {
        distance[i] =  -1;
    }

    //BFS queue for the pairs
    int qs = 0;
    int qe = 0;

    for (int i = 0; i < N; ++i) {
        int id = Id(i, i);
        distance[id] = 0;
        que[qe++] = id;
    }

    //there are more nodes in the queue
    while (qs < qe && qe < noOfPair) {
        int q_id = que[qs++];
        int q_dist = distance[q_id];

        //will process the pair with id q_id now
        int q_s1 = s1fromId(q_id); //the first state in the pair
        int q_s2 = s2fromId(q_id, q_s1); //the second state in the pair (we are sure that q_s1 >= q_s2)


        for (int p = 0; p < P; p++) {
            int* p_ia = &ia[p * N]; //this is the inverse automata for letter p
            int* p_iap = &iap[p * (N + 1)]; //and its state pointers
            int iap_s1_limit = p_iap[q_s1 + 1];
            int iap_s2_limit = p_iap[q_s2 + 1];

            if (p_iap[q_s2] == iap_s2_limit) continue;

            for (int iap_s1_ptr = p_iap[q_s1]; iap_s1_ptr < iap_s1_limit; ++iap_s1_ptr) {
                int ia_s1 = p_ia[iap_s1_ptr];
                for (int iap_s2_ptr = p_iap[q_s2]; iap_s2_ptr < iap_s2_limit; ++iap_s2_ptr) {
                    int ia_s2 = p_ia[iap_s2_ptr];
                    int ia_id = Id(ia_s1, ia_s2);

                    if (distance[ia_id] < 0) { //we found an unvisited pair. so we need to add this to the queue
                        distance[ia_id] = q_dist + 1;
                        next[ia_id] = q_id;
                        letter[ia_id] = p;
                        que[qe++] = ia_id;
                    }
                }
            }
        }
    }

    int* levels = new int[200];
    memset(levels, 0, 200 * sizeof(int));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            int id = Id(i, j);

            if (distance[id] == -1) {
                for (int p = 0; p < P; p++) {
                    int ts1 = a[p * N + i];
                    int ts2 = a[p * N + j];

                    int tid = Id(ts1, ts2);
                    //out << "tid " <<  tid << ": distance is " << distance[tid] << endl;
                }

                //cout << "automata is not synchronizing. pair " << id << " - (" << i << ", " << j << ") is not mergable\n";
                delete [] levels;
                delete [] distance;
                delete [] letter;
                delete [] que;
                delete [] next;
                return false;
            } else {
                levels[distance[id]]++;
            }
        }
    }

    delete [] levels;
    delete [] distance;
    delete [] letter;
    delete [] que;
    delete [] next;
    return true;
}

void printAutomata(int* automata, int N, int p) {
    cout << "Automata ----------------------" << endl;
    for (int i = 0; i < p; ++i) {
        cout << "letter " << (char)(i + 97) << ":\t";
        for (int j = 0; j < N; ++j) {
            cout << automata[i * N + j] << "\t";
        }
        cout << endl;
    }
}

void writeAutomata(int *automata,int N,int P,fstream &out){
    out << "Automata ----------------------" << endl;
    for (int i = 0; i < P; ++i) {
        out << "letter " << (char)(i + 97) << ":\t";
        for (int j = 0; j < N; ++j) {
            out << automata[i * N + j] << "\t";
        }
        out << endl;
    }
}

unsigned long long int shortestPath(int* a, int N, int P, PNode* &path) {
    unsigned long long int noOfNodes = pow (2.0, double(N));
    unsigned long long int* distance = new unsigned long long int[noOfNodes];
    unsigned long long int* prev = new unsigned long long int[noOfNodes];
    unsigned long long int* letter = new unsigned long long int[noOfNodes];
    priority_queue< DistID<unsigned long long int> > que; //ege: min heap for the weighted tree

    for (unsigned long long int i = 0; i < noOfNodes; i++) {
        distance[i] =  ULLONG_MAX;
    }

    //dijkstra queue for the pairs
    distance[noOfNodes-1] = 0;
    prev[noOfNodes - 1] = noOfNodes - 1;
    que.push(DistID<unsigned long long int>(noOfNodes - 1, 0, 0));

    unsigned long long int * q_sN = new unsigned long long int[N];
    unsigned long long int * nextState = new unsigned long long int[N];
    //there are more nodes in the queue
    unsigned long long int counter = 1;
    while (!que.empty()) {
        unsigned long long int q_id = (que.top()).id;
        if (q_id < 0) cout << q_id << endl;
        que.pop();
        unsigned long long int q_dist = distance[q_id];
        unsigned long long int temp_q_id = q_id;
        //int bin_id = 0;

        //will process the pair with id q_id now
        for(unsigned long long int i = 0; i < N ; i++)
        {
            q_sN[i] = temp_q_id % 2;
            temp_q_id = temp_q_id >> 1;
            //if(q_sN[i] != 0)
            //	bin_id += pow(10,N-i-1);
        }
#ifdef DEBUG
        //cout << "will process\t" << bin_id << "\twith id\t" << q_id << "\twith distance\t" << q_dist << endl;
#endif
        for (unsigned long long int p = 0; p < P; p++) {
            //unsigned long long int p = orderedInputs[j];
            memset(nextState, 0, sizeof(unsigned long long int) * N);

            for(unsigned long long int i = 0; i < N; i++)
            {
                if(q_sN[i] != 0)
                {
                    nextState[a[p*N+i]] = 1;
                }

            }
            unsigned long long int id = 0;
            for(unsigned long long int i = 0; i < N; i++)
            {
                if(nextState[i] == 1)
                {
                    id += (1 << i);
                }
            }

            if (distance[id] > q_dist + 1) //berk: key change
            {
                distance[id] = q_dist + 1;
                prev[id] = q_id;
                letter[id] = p;
                que.push(DistID<unsigned long long int>(id,distance[id], counter++));
            }
        }
    }

    delete[] q_sN;
    delete[] nextState;

    unsigned long long int mindist = ULLONG_MAX;
    unsigned long long int minid;
    for (unsigned long long int i = 0; i < N; i++) {
        unsigned long long int pw = pow(2, i);
        if (distance[pw] < mindist) {
            mindist = distance[pw];
            minid = pw;
        }
    }

    unsigned long long int length = 0;
    unsigned long long int s = minid;
    unsigned long long int limit = pow(2, N) - 1;
    vector<char> seq;
    while (s < limit) {
        unsigned long long int p = letter[s];
        seq.push_back(char(p + 97));
        s = prev[s];
        length++;
    }

    // cout << "Shortest Path: ";
    // for (long long int i = seq.size() - 1; i >= 0; i--) {
    //     cout << seq[i] << " ";
    // }
    //cout << endl << "Shortest Path Length:" << length << endl;
    delete [] distance;
    delete [] prev;
    delete [] letter;
    return length;
}

bool checkInverse(int *a, int* iap, int* ia, int N, int P) {
    for (int p = 0; p < P; p++) {
        for (int i = 0; i < N; i++) {
            int target = a[p * N + i];

            int found = 0;
            for (int iaptr = iap[p * (N + 1) + target]; iaptr < iap[p * (N + 1) + target + 1]; ++iaptr) {
                int incoming = ia[p * N + iaptr];
                if (i == incoming) {
                    found = 1;
                    break;
                }
            }

            if (!found) {
                cout << "something is wrong " << i << " goes to " << target << " with " << p << " but it is not in the inverse automata\n";
                return false;
            }
        }
    }

    for (int p = 0; p < P; p++) {
        for (int i = 0; i < N; i++) {
            for (int iaptr = iap[p * (N + 1) + i]; iaptr < iap[p * (N + 1) + i + 1]; ++iaptr) {
                int source = ia[p * N + iaptr];
                if (a[p * N + source] != i) {
                    cout << "something is wrong " << i << " has " << source << " in inverse automata but it " << source << " goes to " << a[p * N + source] << " with " << p << "\n";
                    return false;
                }
            }
        }
    }

    return true;
}

struct Auto{
    int * automata= nullptr;
    unsigned long long int  score=0;
    long long int stats[7] = {0,0,0,0,0,0,0};
    /*
    ~Auto(){
        delete[] automata;
    }
     */
};

bool comparison(Auto a, Auto b){
    return (a.score >  b.score);
}


// BELOW ARE DIFFERENT MUTATION TECHNIQUES ------------------------------------------------------

void mutateSelfLoop(int *ptr,int loc, int N){
    ptr[loc] = (loc%N);
};

void mutateInDegree(int *ptr, int N, int P,int j){
    vector<int> indegs(N);
    int opt = 0;
    if( j >= N) { opt = 1;} // meaning its B transition
    else if(j<N) { opt = 2;} // meaning its A transitions

    // go for A's or B's
    for(int i = N*opt - N;i<N*opt;i++){
        indegs[ptr[i]]++;
    }

    vector<double> probs(N);
    for(int i=0; i <  probs.size();i++){
        if(indegs[i] == 0){ probs[i] = 1.0;}
        else{
            probs[i] = 1/(double)indegs[i];
        }
    }
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::sort(probs.begin(),probs.end());
    double trial = dis(gen);
    for(int i = 0;i<probs.size(); i++){
        if(trial< probs[i]){
            ptr[j] = i;
            break;
        }
    }
}
// ------------------------------------------------------------------------------------------------

// BELOW ARE DIFFERENT CROSS OVER TECHNIQUES -----------------------------------------------------
void regular_threeway(int N,int  P,int K,int a_amount,vector<Auto> &alist){
    random_device ad;
    mt19937 generator(ad());
    std::uniform_int_distribution<> dis4(0,a_amount-1);
    std::uniform_int_distribution<> special(0,(N*P)-1);

    int f,s,t;
    f = dis4(generator);
    s = dis4(generator);
    t = dis4(generator);

    while(f == s && f == t){
        f = dis4(generator);
        s = dis4(generator);
        t = dis4(generator);
    }

    int * fir_parent = alist[f].automata;
    int * sec_parent = alist[s].automata;
    int *tir_parent = alist[t].automata;
    int i_loc;
    i_loc = special(generator);
    if(i_loc == (N*P)-1){ i_loc--;}
    int ii_loc = special(generator);
    while(ii_loc+i_loc > N*P){
        ii_loc = special(generator);
    }
    int iii_loc = special(generator);
    while(ii_loc+i_loc+iii_loc > N*P){
        iii_loc = special(generator);
    }

    for(int i = 0; i<i_loc;i++){
        alist[K].automata[i] = fir_parent[i];
    }
    for(int i=i_loc; i<i_loc+ii_loc ; i++){
        alist[K].automata[i] = sec_parent[i];
    }
    for(int i=i_loc+ii_loc ; i< N*P;i++){
        alist[K].automata[i] = tir_parent[i];
    }

}

void node_wise_mix(int N,int P,int K,int a_amount,vector<Auto> &alist){
    random_device ad;
    mt19937 generator(ad());
    std::uniform_int_distribution<> dis3(0,N-1);
    std::uniform_int_distribution<> dis4(0,a_amount-1);
    int *target = alist[K].automata;
    for(int i = 0; i<N; i ++){
        int pop_tar = dis4(generator);
        int *poptar = alist[pop_tar].automata;
        int node = dis3(generator);
        target[i] = poptar[node];
        target[i+N] = poptar[node+N];
    }

}

void alphabet_wise_mix(int N,int P,int K,int a_amount,vector<Auto> &alist){
    random_device ad;
    mt19937 generator(ad());
    std::uniform_int_distribution<> dis4(0,a_amount-1);
    int *target = alist[K].automata;
    int par1 = dis4(generator);
    int par2 = dis4(generator);

    // ASSUMING P IS 2
    for(int i = 0; i<N;i++){
        target[i] = alist[par1].automata[i];
    }
    for(int i = N; i<N*P;i++){
        target[i] = alist[par2].automata[i];
    }

}

void complete_random(int N,int P,int K,int a_amount,vector<Auto> &alist){
    random_device ad;
    mt19937 gen(ad());
    std::uniform_int_distribution<> dis(0,N-1);
    int *target = alist[K].automata;
    for (int i = 0; i < P * N; ++i) {
        target[i] =  dis(gen);
    }

}


void create_cerny(int N, int P , int *&target){

    for(int i = 0 ; i < N ; i++){
        target[i] = (i+1) % N;

    }
    for(int i=N; i<N*P;i++){
        if(i == (N*P-1)){
            target[i] = 0;
            break;
        }
        target[i] = i%N;
    }

}


// -----------------------------------------------------------------------------------------------


int main(int argc,  char** argv)
{
    int N,P;
    unsigned int seed;
    int auto_amounts;
    unsigned long long int cycles;
    int ratio;

    /*
    if (argc != 6)
    {
        cout << "Usage: " << argv[0] << " no_states, alphabet_size, population size,cylces, extinction ratio \n" << endl;
        return 1;
    }
    N = atoi(argv[1]);
    P = atoi(argv[2]);
    auto_amounts = atoi(argv[3]);
    cycles = strtoull(argv[4],nullptr,10);
    ratio = atoi(argv[5]);

*/
    cin >> N >> P >> auto_amounts >> cycles >> ratio;
    int high_ = (N-1)*(N-1);
    int low = high_ - floor(high_*45/100);
    int NUM_OF_ELITES = 10;

    int *inv_automata_ptrs = new int[P * (N + 1)];
    int *inv_automata = new int[P * N];

    vector<Auto> automata_list;

    string file_name = "genetic_log_" + to_string(N) + "_" + to_string(P) +".txt";
    string file_name2  = "high_points_" + to_string(N) + "_" + to_string(P)  + ".txt";
    string file_name3  = "elites" + to_string(N) + "_" + to_string(P)  + ".txt";
    fstream out(file_name3, fstream::in | fstream::out | fstream::app);
    if(!out){
        fstream out(file_name3, fstream::in | fstream::out | fstream::trunc);
    }


    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0,N-1);
    std::uniform_int_distribution<> dis2(0,8);
    std::uniform_int_distribution<> dis3(0,5);



    while(automata_list.size() < auto_amounts) {
        int * automata = new int[N*P];

        for (int i = 0; i < P * N; ++i) {
            automata[i] =  dis(gen);
        }

        // black magic below
        /*
        int* inv_automata_ptrs = new int[P * (N + 1)];
        int* inv_automata = new int[P * N];
        */
        for (int i = 0; i < P; ++i) {
            int *a = &(automata[i * N]);
            int* ia = &(inv_automata[i * N]);
            int* iap = &(inv_automata_ptrs[i * (N + 1)]);

            memset(iap, 0, sizeof(int) * (N + 1));
            for (int j = 0; j < N; j++) {iap[a[j] + 1]++;}
            for (int j = 1; j <= N; j++) {iap[j] += iap[j - 1];}
            for (int j = 0; j < N; j++) {ia[iap[a[j]]++] = j;}
            for (int j = N; j > 0; j--) {iap[j] = iap[j - 1];} iap[0] = 0;
        }

        // check if these things return true
        // TODO: Error below
        bool isSync = syncCheck(automata, inv_automata_ptrs, inv_automata, N, P);
        // if it fits our criteria, put it into the list
        if(isSync){
            PNode *path = nullptr;
            unsigned long long int shortest = shortestPath(automata, N, P, path);

            Auto temp;
            temp.automata = automata;
            temp.score = shortest;
            automata_list.push_back(temp);
            //cout << "One added to the list... \n";
        } else {
            delete[] automata;
        }

#ifdef DEBUG
        // printAutomata(automata, N, P);
#endif

        // it is not synchronizing and/or doesnt have inverse

        /*
        delete[] inv_automata_ptrs;
        delete[] inv_automata;
        */
    }
    /*
    cout << "The population for the genetic algorithm has been created... \n";

    out << "Randomly created population:\n";
    out << "With parameters Nodes: " << N << " Alphabet: " << P<< " Population Size: "<< auto_amounts<<" Cut ratio: 1/"<< ratio <<"\n";
    unsigned long long int ftotal_score = 0;
    for (int i =0 ; i< automata_list.size();i++){
        out << automata_list[i].score << " ";
        ftotal_score += automata_list[i].score;
    }

    out << "Average: " << ftotal_score/automata_list.size() << "\n";
    */

    // begin genetic cycle
    // can easily migrate this into a function
    unsigned long long int begining = 0;
    unsigned long long int Generation_ = 1;


    // throw in a cerny
    int *ptr = automata_list[0].automata;
    //create_cerny(N,P,ptr);
    PNode *ppath = nullptr;
    unsigned long long int shortest = shortestPath(ptr, N, P, ppath);
    ppath = nullptr;
    automata_list[0].score = shortest;



    while(cycles){

        if(Generation_ % 1000 == 0){
            cout << "Currently in Generation: " << Generation_<<endl;
            // out<< "Generation: " << Generation_ << endl;
            /*
            for (int i = 0;i < NUM_OF_ELITES;i++){
                out << automata_list[i].score << " ";
            }
            out << "\n";
             */
        }
        // sort the list


        if(Generation_ == 1){
            begining = automata_list[0].score;
        }

        int modif = (int)floor((auto_amounts*(ratio-1))/ratio);



        // MUTATION LOGIC BELOW
        // mutate < ------- > mutate

        for(int i=0;i < automata_list.size();) {
            vector<int> re(N*P);

            if (dis3(gen) == 5 || dis3(gen) ==4 || dis3(gen) == 1) {
                // will mutate now
                int indic;
                int *ptr = automata_list[i].automata;
                for(int j = 0 ; j< N*P;j++){
                    int r = dis2(gen);

                    // ADD THE MUTATION LOGICS BELOW WITH THEIR PROBABILITIES **
                    if(r == 5 || r== 7 || r==0){
                        re[j] = ptr[j];
                        ptr[j] = dis(gen);
                        indic = RAND_MUT;
                    }
                    else if(r==1 || r==4 ){
                        re[j] = ptr[j];
                        // mutation method
                        //mutateSelfLoop(ptr,j,P);
                        mutateInDegree(ptr,N,P,j);
                        indic = MUTATE_INDEG;
                    }

                    else if(r==2 ){
                        re[j] = ptr[j];
                        mutateSelfLoop(ptr,j,P);
                        indic = MUTATE_SELF;
                    }

                    // MUTATION LOGICS ABOVE *************************************

                }
                /*
                int *inv_automata_ptrs = new int[P * (N + 1)];
                int *inv_automata = new int[P * N];
                */
                for (int k = 0; k < P; ++k) {
                    int *a = &(ptr[k * N]);
                    int* ia = &(inv_automata[k * N]);
                    int* iap = &(inv_automata_ptrs[k * (N + 1)]);

                    memset(iap, 0, sizeof(int) * (N + 1));
                    for (int z = 0; z < N; z++) {iap[a[z] + 1]++;}
                    for (int z = 1; z <= N; z++) {iap[z] += iap[z - 1];}
                    for (int z = 0; z < N; z++) {ia[iap[a[z]]++] = z;}
                    for (int z = N; z > 0; z--) {iap[z] = iap[z - 1];} iap[0] = 0;
                }
                PNode *path = nullptr;
                // TODO: CHECK IF SYNC ...
                bool isSync = syncCheck(ptr, inv_automata_ptrs, inv_automata, N, P);
                if(isSync) {
                    unsigned long long int shortest = shortestPath(automata_list[i].automata, N, P, path);
                    automata_list[i].score = shortest;
                    ptr = nullptr;
                    automata_list[i].stats[indic]++;
                    i++;
                }else{
                    for(int z=0;z<N*P;z++){
                        // revert back
                        if(re[z] != 0) ptr[z] = re[z];
                    }
                }
            }
            else{i++;}

        }
        sort(automata_list.begin(),automata_list.end(),comparison);

        // cross < ---  > over
        for(int K=modif; K < automata_list.size();){

            // use the cross over technique you want
            int indic_cross;
            // regular_threeway( N,P,K,auto_amounts,automata_list);
            int xx = dis2(gen);
            // CROSS-OVER LOGIC IN THE FOLLOWING BOX ****************************************
            if(xx == 1 ) {complete_random(N,P,K,auto_amounts,automata_list);indic_cross=COMP_RAND;}
            else if (xx == 2 || xx == 3) {regular_threeway(N,P,K,auto_amounts,automata_list);indic_cross=THREE_CROSS;}
            else if(xx==5 || xx==4 ){node_wise_mix(N,P,K,auto_amounts,automata_list); indic_cross=NODE_WISE;}
            else {alphabet_wise_mix(N,P,K,auto_amounts,automata_list);indic_cross = ALPH_WISE;}
            // CROSS-OVER LOGIC ABOVE *************indic_cross******************************************



            // black magic below
            /*
            int *inv_automata_ptrs = new int[P * (N + 1)];
            int *inv_automata = new int[P * N];
            */
            for (int i = 0; i < P; ++i) {
                int *a = &(automata_list[K].automata[i * N]);
                int *ia = &(inv_automata[i * N]);
                int *iap = &(inv_automata_ptrs[i * (N + 1)]);

                memset(iap, 0, sizeof(int) * (N + 1));
                for (int j = 0; j < N; j++) { iap[a[j] + 1]++; }
                for (int j = 1; j <= N; j++) { iap[j] += iap[j - 1]; }
                for (int j = 0; j < N; j++) { ia[iap[a[j]]++] = j; }
                for (int j = N; j > 0; j--) { iap[j] = iap[j - 1]; }
                iap[0] = 0;
            }

            // check if these things return true
            bool isSync = syncCheck(automata_list[K].automata, inv_automata_ptrs, inv_automata, N, P);
            // if it fits our criteria, put it into the list
            if(isSync){
                PNode *path = nullptr;
                unsigned long long int shortest = shortestPath(automata_list[K].automata, N, P, path);
                automata_list[K].score = shortest;
                for(int x= 0; x < 7 ; x++){
                    automata_list[K].stats[x] = 0;
                }
                automata_list[K].stats[indic_cross]++;
                K++;

            }
        }
        cycles--;
        Generation_++;

        /*
       if(Generation_ %5000 == 0) {
            out << " Generation: " << Generation_ << "\n";
            unsigned long long int total_score = 0;

            for (int i = 0; i < automata_list.size(); i++) {
                out << automata_list[i].score << " ";
                total_score += automata_list[i].score;
            }
            out << " Average: " << total_score / automata_list.size() << "\n";

            out << " The automatas: \n";
            for (int i = 0; i < automata_list.size(); i++) {
                out << "Automata number " << i + 1 << ":\n";
                writeAutomata(automata_list[i].automata, N, P, out);
            }


            out << "============================================================================ \n";
        }
*/

        /*
        for(int ii=0; ii< NUM_OF_ELITES; ii++){
            if(automata_list[ii].score > low){
                int rndm = dis2(gen);

                // take with probabilty from elites to let the growth
                if(rndm == 5){
                    writeAutomata(automata_list[ii].automata,N,P,out);
                    out  << "score: " << automata_list[ii].score << "\n";

                    bool check = true;

                    while(check) {
                        // clean the automata
                        complete_random(N, P, ii, auto_amounts, automata_list);

                        //int *inv_automata_ptrs = new int[P * (N + 1)];
                        //int *inv_automata = new int[P * N];

                        for (int i = 0; i < P; ++i) {
                            int *a = &(automata_list[ii].automata[i * N]);
                            int *ia = &(inv_automata[i * N]);
                            int *iap = &(inv_automata_ptrs[i * (N + 1)]);

                            memset(iap, 0, sizeof(int) * (N + 1));
                            for (int j = 0; j < N; j++) { iap[a[j] + 1]++; }
                            for (int j = 1; j <= N; j++) { iap[j] += iap[j - 1]; }
                            for (int j = 0; j < N; j++) { ia[iap[a[j]]++] = j; }
                            for (int j = N; j > 0; j--) { iap[j] = iap[j - 1]; }
                            iap[0] = 0;
                        }

                        // check if these things return true
                        bool isSync = syncCheck(automata_list[ii].automata, inv_automata_ptrs, inv_automata, N, P);
                        // if it fits our criteria, put it into the list
                        if(isSync){
                            PNode *path = nullptr;
                            unsigned long long int shortest = shortestPath(automata_list[ii].automata, N, P, path);
                            automata_list[ii].score = shortest;
                            check = false;

                        }

                        //delete[] inv_automata;
                        //delete[] inv_automata_ptrs;

                         }

                }
            }
        }
        */



    }
    /*
    sort(automata_list.begin(),automata_list.end(),comparison);

    unsigned long long int final = automata_list[0].score;
    out << begining << " " << final << "\n";
*/
    /*#define MUTATE_SELF 0;
#define MUTATE_INDEG 1;
#define THREE_CROSS 2;
#define NODE_WISE 3;
#define COMP_RAND 4;
#define RAND_MUT 5;
*/
    std::string dict[7] ={"MUTATE_SELF","MUTATE_INDEG","THREE_CROSS","NODE_WISE","COMP_RAND","RAND_MUT","ALPH_WISE"};
    // print what happened to top 10
    for(int i = 0; i< 10 ;i++){
        out << automata_list[i].score << "\n";
        for(int j = 0; j < 7;j++){
            out  <<dict[j] << ": " << automata_list[i].stats[j] << ", ";
        }
        out << "\n";
        out << " ============================================== \n";
    }
    /*
     *  Delete everything below here
     */
    for(int i = 0 ; i<automata_list.size();i++){
        delete [] automata_list[i].automata;
    }
    delete[] inv_automata;
    delete[] inv_automata_ptrs;
    return 0;
}

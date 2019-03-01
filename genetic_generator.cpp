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

// #define DEBUG

#define Id(s1, s2) ((s1 > s2)?(((s1 * (s1 + 1))/2) + s2):(((s2 * (s2 + 1))/2) + s1)) //this is how we compute the ids
#define IdNaive(s1, s2, N) ((s1 > s2)?((s1 * N) + s2):((s2 * N) + s1)) //this is how we compute the ids
#define s1fromId(id) ((int)(sqrt((2.0 * id) +1.0) - 0.5));
#define s2fromId(id, s1) (id - ((s1 * (s1 + 1))/2));

using namespace std;

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
                return false;
            } else {
                levels[distance[id]]++;
            }
        }
    }
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

void writeAutomata(int *automata,int N,int P,ofstream &out){
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
    /*
    ~Auto(){
        delete[] automata;
    }
     */
};

bool comparison(Auto a, Auto b){
    return (a.score >  b.score);
}
// TODO: HAS INVERSE CHECKE GEREK YOK
// TODO: Add logging procedure to observe the changes to the population


void mutateSelfLoop(int *ptr,int loc, int P){
    ptr[loc] = int(floor(loc/P));
};

void mutateInDegree(int *ptr, int N, int P,int j){
    vector<int> indegs(N);

    // calculate the indegs for each Node
    for(int i=0;i<N*P;i++){
        indegs[ptr[i]] +=1;
    }
    // sort(indegs.begin(),indegs.end());

    random_device ran;
    mt19937 gn(ran());
    std::uniform_int_distribution<> diss(0,N*P);

    int min = N*P*2;
    int loc;
    for(int i = 0;i<indegs.size();i++){
        if(indegs[i]<=min) {
            min = indegs[i];
            loc = i;
        }
    }
    ptr[j] = loc;
}

int main(int argc,  char** argv)
{
    int N,P;
    unsigned int seed;
    int auto_amounts;
    int cycles;
    int ratio;
    
    if (argc != 6)
    {
        cout << "Usage: " << argv[0] << " no_states, alphabet_size, population size,cylces,extinction ratio \n" << endl;
        return 1;
    }
    N = atoi(argv[1]);
    P = atoi(argv[2]);
    auto_amounts = atoi(argv[3]);
    cycles = atoi(argv[4]);
    ratio = atoi(argv[5]);


   //  cin >> N >> P >> auto_amounts >> cycles >> ratio;



    vector<Auto> automata_list;

    string file_name = "genetic_log_" + to_string(N) + "_" + to_string(P) +".txt";
    fstream out(file_name, fstream::in | fstream::out | fstream::app);
    if(!out){
        fstream out(file_name, fstream::in | fstream::out | fstream::trunc);
    }

    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0,N-1);
    std::uniform_int_distribution<> dis2(0,8);
    std::uniform_int_distribution<> dis3(0,4);
    std::uniform_int_distribution<> special(0,(N*P)-1);
    std::uniform_int_distribution<> dis4(0,auto_amounts-1);


    while(automata_list.size() < auto_amounts) {
        int * automata = new int[N*P];


        for (int i = 0; i < P * N; ++i) {
            automata[i] =  dis(gen);
        }




        // black magic below

        int* inv_automata_ptrs = new int[P * (N + 1)];
        int* inv_automata = new int[P * N];

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
           // cout << "One added to the list... \n";
        }
        else{
            // it is not synchronizing and/or doesnt have inverse
            delete[] automata;
            delete[] inv_automata_ptrs;
            delete[] inv_automata;

        }



    }

   // cout << "The population for the genetic algorithm has been created... \n";
    /*
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
    int Generation_ = 1;
    while(cycles){
     //   cout << "Currently in Generation: " << Generation_<<endl;
        // sort the list

        if(Generation_ == 1){
            begining = automata_list[0].score;
        }

        int modif = (int)floor((auto_amounts*(ratio-1))/ratio);

        // eliminate the unwanted population
        //automata_list.erase(automata_list.begin()+modif,automata_list.end());


        // MUTATION LOGIC BELOW
        // mutate < ------- > mutate
        // TODO : Below goes into infinite so put checks...
        for(int i=1;i < automata_list.size();) {
            vector<int> re(N*P);
            if (dis3(gen) == 4) {
                // will mutate now
                int *ptr = automata_list[i].automata;
                for(int j = 0 ; j< N*P;j++){
                    int r = dis2(gen);
                  //  if(r == 5 || r== 7 || r==0 || r==8||r==3){
                       re[j] = ptr[j];
                        mutateInDegree(ptr,N,P,j);
                       // ptr[j] = dis(gen);
                  //  }
                    /*
                        else if(r==1 || r ==2 || r==6||r==4){
                        re[j] = ptr[j];
                        //mutateSelfLoop(ptr,j,P);
                       
                    }
                    */
                    /*
                    else if(r==3 || r==4 ){
                        re[j] = ptr[j];
                        mutateInDegree(ptr,N,P,j);
                    }
                    */
                }
                int *inv_automata_ptrs = new int[P * (N + 1)];
                int *inv_automata = new int[P * N];

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
            int fir,sec;
            fir = dis4(gen);
            sec = dis4(gen);
            if (fir != sec){

                int * fir_parent = automata_list[fir].automata;
                int * sec_parent = automata_list[sec].automata;

                int i_loc;
                i_loc = (special(gen));
                if(i_loc == (N*P)-1){ i_loc--;}
                

                for(int i = 0; i<i_loc;i++){
                    automata_list[K].automata[i] = fir_parent[i];
                }
                for(int i=i_loc; i<N*P ; i++){
                    automata_list[K].automata[i] = sec_parent[i];
                }


                // black magic below
                int *inv_automata_ptrs = new int[P * (N + 1)];
                int *inv_automata = new int[P * N];

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
                    K++;

                }
                else{
                    // it is not synchronizing

                    delete[] inv_automata;
                    delete[] inv_automata_ptrs;
                }

            }
        }
        cycles--;
        Generation_++;
        /*
        out << " Generation: " << Generation_ << "\n";
        unsigned long long int total_score = 0;

        for (int i=0;i < automata_list.size();i++){
            out << automata_list[i].score << " " ;
            total_score += automata_list[i].score;
        }
        out << " Average: " << total_score/automata_list.size() << "\n";

        out << " The automatas: \n";
        for(int i=0; i< automata_list.size();i++){
            out << "Automata number " << i+1 << ":\n";
            writeAutomata(automata_list[i].automata,N,P,out);
        }


        out <<"============================================================================ \n";
         */
    }
    sort(automata_list.begin(),automata_list.end(),comparison);

    unsigned long long int final = automata_list[0].score;
    out << begining << " " << final << "\n";

    /*
     *  Delete everything below here
     */
    for(int i = 0 ; i<automata_list.size();i++){
        delete [] automata_list[i].automata;
    }
    return 0;
}

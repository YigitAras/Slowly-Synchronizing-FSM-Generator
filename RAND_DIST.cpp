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


int main(int argc, char** argv) {
    int N;
    int P = 2; // a and b


    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " no_states, alphabet_size, population size,cylces,extinction ratio \n" << endl;
        return 1;
    }
    N = atoi(argv[1]);
    P = atoi(argv[2]);



    int *inv_automata_ptrs = new int[P * (N + 1)];
    int *inv_automata = new int[P * N];
    string file_name =  "random_dist_" + to_string(N) + "_" + to_string(P)+".txt";

    fstream out(file_name,fstream::in | fstream::out | fstream::app);

    if(!out){
        fstream out(file_name, fstream::in | fstream::out | fstream::trunc);
    }

    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0,N-1);

    int *ex_machina  = new int[N*P];
    memset(ex_machina,0,sizeof(int)*N*P);

    long long int CTR_AMNT = 50000000;
    while(CTR_AMNT){
        // randomly generate the automata
        for (int i = 0; i < P * N; ++i) {
            ex_machina[i] =  dis(gen);
        }
        for (int i = 0; i < P; ++i) {
            int *a = &(ex_machina[i * N]);
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
        bool isSync = syncCheck(ex_machina, inv_automata_ptrs, inv_automata, N, P);
        // if it fits our criteria, put it into the list
        if(isSync){
            PNode *path = nullptr;
            unsigned long long int shortest = shortestPath(ex_machina, N, P, path);
            CTR_AMNT--;
            out << shortest << " ";
        }


    }

    return 0;
}

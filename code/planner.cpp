#include <math.h>
#include <vector>
#include <iostream>
#include <memory>
#include <queue>
#include <utility>
#include <time.h>
#include <unordered_set>
#include <unordered_map>
#include <mex.h>
#include <chrono>
#include <cstdio>
#include <fstream>

using namespace std;
/* Input Arguments */
#define    MAP_IN                  prhs[0]
#define    ROBOT_IN                prhs[1]
#define    TARGET_TRAJ             prhs[2]
#define    TARGET_POS              prhs[3]
#define    CURR_TIME               prhs[4]
#define    COLLISION_THRESH        prhs[5]

/* Output Arguments */
#define    ACTION_OUT              plhs[0]

//access to the map is shifted to account for 0-based indexing in the map, whereas
//1-based indexing in matlab (so, robotpose and goalpose are 1-indexed)
#define GETMAPINDEX(X, Y, XSIZE, YSIZE) ((Y-1)*XSIZE + (X-1))

#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

#define NUMOFDIRS 9
#define NUMOFDIRS_H 8
#define MAX_VAL __DBL_MAX__

ostream& operator<<(ostream& os, const vector<int>& vec)
{
    for(auto iter:vec){
        os<<iter<<" ";
    }
    cout<<endl;
    return os;
}

// initialisation of some global variables which will be used by multiple functions
int x_len, y_len, total_len, total_time, current = 0;
bool created = false;
vector<pair<int,int> > trajectory;

// struct which represents the information of each state in the 3D A* search
struct Node{
    public:
        double g,h;
        bool visited;
        vector<int> parent;
        Node() : g(MAX_VAL), h(-1), parent({-1,-1}), visited(false){
            
        }
};

// struct used to create arrays needed to perform forward 2D A* search
struct forward_heuristic_map{
    vector<bool> visited;
    vector<int> count;
    vector<double> cost;
    forward_heuristic_map(){
        visited.resize((x_len+1)*(y_len+1), false);
        count.resize((x_len+1)*(y_len+1), INT32_MAX);
        cost.resize((x_len+1)*(y_len+1), __DBL_MAX__);
    }
};

// function to convert 2D coordinates to index of 1D array
int get_index(int x, int y){
    return (y-1)*x_len+(x-1);
}


// custom comparator for the open list
class Comparison{
    public:
    bool operator () (pair<double,vector<int> >& a, pair<double, vector<int> >& b){
        return a.first > b.first;
    }
};

// function to check if a move is valid in the grid
bool validate(int curr_x, int curr_y, int dx, int dy){
    if((curr_x+dx >0 && curr_x+dx <= x_len) && (curr_y+dy >0 && curr_y+dy <= y_len))
        return true;
    return false;
}


// class which implements the hash function for unordered map with key as <x,y,t>
class Hash{
    public:
      // as hash function.
        size_t operator()(const vector<int>& p) const
        {
            //y = row x=col t=depth
            return p[1]*x_len*total_time+p[0]*total_time+p[2];
        }
};

// hash function for unordred map with key as <x,y>
class Hash_h{
    public:
      // as hash function.
        size_t operator()(const vector<int>& p) const
        {
            //y = row x=col t=depth
            return p[1]*x_len+p[0];
        }
};

// forward 2D A* search to reach all target states
void forward_A_star(double* map,
    int collision_thresh,
    int x_size,
    int y_size,
    int robotposeX,
    int robotposeY,
    int target_steps,
    double* target_traj,
    vector<bool>& visited,
    vector<int>& count,
    vector<double>& cost,
    bool flag
    )
{
    // initialise data structures
    x_len = x_size;
    y_len = y_size;
    total_len = x_len*y_len;
    total_time = target_steps;
    priority_queue<pair<double,vector<int> >, vector<pair<double, vector<int> > >, Comparison> open_list;
    unordered_set<vector<int>, Hash_h> terminal_states;
    int dX[NUMOFDIRS_H] = {-1, -1, -1,  0,  0,  1, 1, 1};
    int dY[NUMOFDIRS_H] = {-1,  0,  1, -1,  1, -1, 0, 1};
    vector<int> goal = {robotposeX, robotposeY};

    // start timer
    cout<<"Total length of the array is "<<(y_len+1)*(x_len+1)<<endl;
    auto start = chrono::high_resolution_clock::now();
    
    // push all terminal states into the open list
    for(int i=0; i<target_steps; i++){
        terminal_states.insert({(int)target_traj[i], (int)target_traj[target_steps+i]});
    }

    // push starting node to openlist and initialise there cost and count values
    open_list.push({0, {robotposeX, robotposeY}});
    auto top = open_list.top();
    auto curr = top.second;
    int curr_index = get_index(robotposeX, robotposeY);
    count[curr_index] = 0;
    cost[curr_index] = 0;
    
    while(!terminal_states.empty() && !open_list.empty()){
        top = open_list.top();
        curr = top.second;

        // remove a visited terminal state
        curr_index = get_index(curr[0], curr[1]);
        if(terminal_states.find(curr)!=terminal_states.end()){
            terminal_states.erase(curr);
        }
        open_list.pop();

        // if current node is not visited then check if neighbors count can be updated
        if(!visited[curr_index]){
            
            for(int i=0; i<NUMOFDIRS_H; i++){
                // double check if this is right
                auto neighbor = curr;
                neighbor[0]+=dX[i];
                neighbor[1]+=dY[i];
                int index = get_index(neighbor[0], neighbor[1]);
                // cout<<curr_index<<" "<<index<<endl;
                bool condition = false;
                // condition = visited[index];
                if(validate(curr[0], curr[1], dX[i], dY[i]))
                {
                    if(!visited[index] && map[index]<collision_thresh){
                        // check if neighbor count can be reduced
        //                cout<<"Next "<<neighbor[0]<<" "<<neighbor[1]<<endl;
                        if(count[index] > count[curr_index] + 1){
                            // update the count and also the cost inurred 
                            count[index] = count[curr_index] + 1;
                            cost[index] = cost[curr_index] + map[index];
                            // cout<<"neighbor visited is "<<neighbor;
                            open_list.push({count[index], neighbor});
                        }
                    }
                    // if neighbor is unreachable dont push into the open list and make visited as true
                    else if(map[index]>collision_thresh){
                        visited[index] = true;
                    }   
                }
                
            }
            visited[curr_index]= true;
        }
    }
    
    /*
    ofstream ofs("output_1.txt");
        for(int i=0; i<target_steps; i++){
            // cout<<get_index(target_traj[i], target_traj[target_steps+i])<<endl;
            ofs<<count[get_index(target_traj[i], target_traj[target_steps+i])]<<endl;
        }
    ofs.close();
    */
    
    auto now = chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
    cout<<"Time taken to compute forward heuristics is : "<<duration.count()<<endl;
    return;
}

// perform reverse A* with terminal states as the initial states to estimate an informed heuristic
void reverse_A_star(
    double* map,
    int collision_thresh,
    int x_size,
    int y_size,
    int robotposeX,
    int robotposeY,
    int target_steps,
    double* target_traj,
    vector<bool>& visited,
    vector<int>& count,
    vector<double>& cost,
    unordered_set<vector<int>,Hash_h>& terminal
    )
{
    // initialise the variables and data structures
    x_len = x_size;
    y_len = y_size;
    total_len = x_len*y_len;
    total_time = target_steps;
    priority_queue<pair<double,vector<int> >, vector<pair<double, vector<int> > >, Comparison> open_list;
    int dX[NUMOFDIRS_H] = {-1, -1, -1,  0,  0,  1, 1, 1};
    int dY[NUMOFDIRS_H] = {-1,  0,  1, -1,  1, -1, 0, 1};
    vector<int> terminal_state = {robotposeX, robotposeY};
    vector<int> goal = {robotposeX, robotposeY};
    int pointer;

    // push all terminal states into the open list
    for(auto iter: terminal){
        // cout<<iter;
        pointer = get_index(iter[0], iter[1]);
        open_list.push({0, iter});
        cost[pointer] = 0;
        count[pointer] = 0;

    }
    // cout<<"Here"<<endl;
    auto top = open_list.top();
    auto curr = top.second;
    int curr_index = get_index(curr[0], curr[1]);
    // cout<<"Searching"<<endl;
    while(!open_list.empty()){
        top = open_list.top();
        curr = top.second;
        open_list.pop();
        curr_index = get_index(curr[0], curr[1]);
        // cout<<"Current"<<endl;
        // cout<<curr<<endl;
        // cout<<curr_index<<endl;
        if(!visited[curr_index]){
            // cout<<curr<<endl;
            // cout<<"neighbors\n";
            for(int i=0; i<NUMOFDIRS_H; i++){
                // double check if this is right
                auto neighbor = curr;
                neighbor[0]+=dX[i];
                neighbor[1]+=dY[i];
                int index = get_index(neighbor[0], neighbor[1]);
                // if valid move and not visited and reachable state
                // cout<<neighbor<<endl;
                // cout<<validate(curr[0], curr[1], dX[i], dY[i])<<endl;
                if(validate(curr[0], curr[1], dX[i], dY[i])){
                    // cout<<index<<endl;
                    // cout<<visited[index]<<endl;
                    // cout<<map[index]<<endl;
                    if(!visited[index] && map[index]<collision_thresh)
                    {
                        // check if neighbor cost can be reduced
                        // cout<<"Next "<<neighbor[0]<<" "<<neighbor[1]<<endl;
                        if(cost[index] > cost[curr_index] + map[index]){
                            // update the cost if it can be reduced
                            // also update the count of steps taken in such a case
                            cost[index] = cost[curr_index] + map[index];
                            count[index] = count[curr_index] + 1;
                            // cout<<"Pushed "<<neighbor<<endl;
                            open_list.push({cost[index], neighbor});
                        }
                    }
                    else if(map[index]>collision_thresh){
                        // cout<<"Entered here"<<endl;
                        visited[index] = true;
                    }
                }
            }
            visited[curr_index] = true;
        }
        // cout<<(terminal.find(curr)==terminal.end() && !open_list.empty())<<endl;
        // if(!(curr!=goal && !open_list.empty()))
        // cout<<(curr!=goal && !open_list.empty())<<endl;
    }
    return;
}

// function to check if we have reached any of the terminal states
bool isTerminal(vector<int> curr, unordered_set<vector<int>, Hash_h>& terminal_states, unordered_map<vector<int>, int, Hash_h>& time_lookup, int thresh){
    if(terminal_states.find({curr[0], curr[1]})!=terminal_states.end() && curr[2]==(time_lookup[{curr[0], curr[1]}]-thresh))
        return true;
    return false;
    
}


// function to select a window of trajectory points towards which we will perform
// reverse A* search and then forward 3D A* search
void select_terminal_states(unordered_set<vector<int>,Hash_h>& terminal, 
    forward_heuristic_map fh_map,
    unordered_map<vector<int>, int, Hash_h>& time_lookup,
    double *target_traj, 
    int target_steps,
    int &thresh,
    int& max_time,
    vector<int>& max_time_state,
    double* map)
    {

    // ofstream logger("terminal.txt");
    // initialise the variables and data structures
    double min_cost_ub = MAX_VAL;
    int index = -1, delta = 0, delta1=0;
    int pointer = 0;
    int window_size = 25;
    // until a minimum upper bound trajectory point can be selected halve the threshold and recheck
    while(true){
        min_cost_ub = MAX_VAL;
        index = -1;
        for(int i=0; i<target_steps; i++){
            pointer = get_index((int)target_traj[i],(int)target_traj[i+target_steps]);
            if(fh_map.count[pointer]>i-thresh)
                continue;
            else{
                // cout<<fh_map.count[pointer]<<" "<<i<<" "<<thresh<<endl;
                delta = i-fh_map.count[pointer];
                double ub = fh_map.cost[pointer]+delta*map[pointer];
                // double lb = Heuristic_Map[temp].g_f;
                double avg = ub;
                if(min_cost_ub > avg){
                    min_cost_ub = avg;
                    index = i;
                }
            }
        }
        // if no such trajectory point threshold as cushion was found, reduce th threshold
        if(index==-1){
            thresh = thresh/2;
        }
        else{
            break;
        }
    }

    vector<int> temp_ = {(int)target_traj[index], (int)target_traj[index+target_steps]};
    // logger<<"The upperbound cost is "<<min_cost_ub<<" "<<index<<" "<<Heuristic_Map[temp_].count<<endl;
    // logger<<"The selected terminal states are "<<endl;
    // cout<<"Terminal states are\n";

    // select a window of trajectory points around the minimum upper bound cost trajectory point identified
    for(int i = max(0,index-window_size); i<=min(target_steps-1,index+window_size); i++){
        vector<int> temp = {(int)target_traj[i], (int)target_traj[i+target_steps]};
        pointer = get_index(temp[0], temp[1]);
        if(fh_map.count[pointer]>i-thresh)
            continue;
        else{
            terminal.insert(temp);
            // cout<<temp<<endl;
            time_lookup[temp] = i;
            max_time = i;
            max_time_state = temp;
        }
    }
}

static void planner(
        double* map,
        int collision_thresh,
        int x_size,
        int y_size,
        int robotposeX,
        int robotposeY,
        int target_steps,
        double* target_traj,
        int targetposeX,
        int targetposeY,
        int curr_time,
        double* action_ptr
){
    cout<<"curr_time is "<<curr_time<<endl;
    cout<<"current is "<<current<<endl;
    cout<<"robot pose is "<<robotposeX<<" "<<robotposeY<<endl;
    cout<<"target pose is "<<targetposeX<<" "<<targetposeY<<endl;
    x_len = x_size;
    y_len = y_size;
    total_len = x_len*y_len;
    total_time = target_steps;
    int thresh = (int)(x_len*y_len)*4e-6;
    cout<<"Threshold "<<thresh<<endl;
    auto start = chrono::high_resolution_clock::now();
    if(created){
        action_ptr[0] = trajectory[current].first;
        action_ptr[1] = trajectory[current].second;
        current++;
        return;
    }
    ofstream points("points.txt");
    static unordered_map<vector<int>,Node, Hash> Graph;
    static priority_queue<pair<double,vector<int> >, vector<pair<double, vector<int> > >, Comparison> open_list;
    int dX[NUMOFDIRS] = {-1, -1, -1,  0,  0,  1, 1, 1, 0};
    int dY[NUMOFDIRS] = {-1,  0,  1, -1,  1, -1, 0, 1, 0};
    static unordered_set<vector<int>,Hash_h> terminal;
    static unordered_map<vector<int>, int, Hash_h> time_lookup;
    static forward_heuristic_map fh_map, rh_map;

    cout<<"Started Forward A* search"<<endl;
    forward_A_star(map, collision_thresh, x_size, y_size, robotposeX, robotposeY, target_steps, target_traj, fh_map.visited, fh_map.count, fh_map.cost, true);
    cout<<"Finished Forward A* search"<<endl;
    auto now = chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
    cout<<"Time taken to compute forward heuristics is : "<<duration.count()<<endl;

    cout<<"Terminal states"<<endl;
    int max_time,count=0;
    vector<int> max_time_state;
    select_terminal_states(terminal, fh_map, time_lookup, target_traj, target_steps, thresh, max_time, max_time_state, map);
    now = chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
    cout<<"Time taken to compute Terminal States is : "<<duration.count()<<endl;
    
    // clear unecessary arrays to save memory
    fh_map.cost.clear();
    fh_map.count.clear();
    fh_map.visited.clear();

    cout<<"Started Reverse A* search"<<endl;
    reverse_A_star(map, collision_thresh, x_size, y_size, robotposeX, robotposeY, target_steps, target_traj, rh_map.visited, rh_map.count, rh_map.cost, terminal);
    cout<<"Finished Reverse A* search"<<endl;
    now = chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
    cout<<"Time taken to compute reverse heuristics is : "<<duration.count()<<endl;
    
    // initialise the openlist with the starting position of the robot
    int curr_index = get_index(robotposeX, robotposeY);
    Graph[{robotposeX, robotposeY, 0}] = Node();
    Graph[{robotposeX,robotposeY, 0}].g = 0;
    Graph[{robotposeX,robotposeY,0}].h = rh_map.cost[curr_index];
    open_list.push({0,{robotposeX, robotposeY, 0}});
    auto top = open_list.top();
    auto curr = top.second;
    
    while(!isTerminal(curr, terminal, time_lookup, thresh) && !open_list.empty()){
        auto now = chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms = now - start;
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);

        if(duration.count()%1000==0){
            if(duration.count()>=(thresh+1)*1000){
                // Terminate if it takes more than threshold seconds
                cout<<"Time out"<<endl;
                return;
            }
            cout<<"Time is "<<duration.count()<<endl;
        }

        top = open_list.top();
        curr = top.second;
        open_list.pop();
        curr_index = get_index(curr[0], curr[1]);

        // identify if expanding the current state is useful or not
        int min_time_req = max(abs(max_time_state[0]-curr[0]), abs(max_time_state[1]-curr[1]));
        int time_diff = max_time-curr[2]-thresh;
        bool condition = time_diff<min_time_req;
        if(Graph[curr].visited||condition){
            continue;
        }
        // make current node visited
        Graph[curr].visited = true;

        // cout<<"Visited f value is "<<Graph[curr].g+Graph[curr].h<<endl;
        // points<<curr[0]<<" "<<curr[1]<<" "<<curr[2]<<" "<<Graph[curr].g<<" "<<Graph[curr].h<<endl;

        // check neighbors
        for(int i=0; i<NUMOFDIRS; i++){
            // double check if this is right
            auto neighbor = curr;
            neighbor[0]+=dX[i];
            neighbor[1]+=dY[i];
            neighbor[2]++;
            int index = get_index(neighbor[0], neighbor[1]);
            // if valid move and not visited and reachable state
            if(validate(curr[0], curr[1], dX[i], dY[i])
                && Graph.find(neighbor) == Graph.end()){
                Graph[neighbor] = Node();
            }
            // cout<<"validate "<<validate(curr[0], curr[1], dX[i], dY[i])<<endl;
            // cout<<"visited "<<Graph[neighbor].visited<<endl;
            // cout<<"reachable "<<( map[index]<collision_thresh)<<endl;
            if(validate(curr[0], curr[1], dX[i], dY[i])
                && !Graph[neighbor].visited
                && map[index]<collision_thresh)
            {
                // check if neighbor cost can be reduced
                // cout<<"Next "<<neighbor[0]<<" "<<neighbor[1]<<" t is: "<<neighbor[2]<<endl;
                // cout<<Graph[neighbor].g<<" "<<Graph[curr].g<<" "<<map[index]<<endl;
                if(Graph[neighbor].g > Graph[curr].g + map[index]){
                    // update the g value
                    // update the f value
                    // update the parent
                    // insert the updated items into the priority queue
                    Graph[neighbor].g = Graph[curr].g + map[index];
                    if(Graph[neighbor].h==-1){
                        if(neighbor[2]<=max_time){
                            Graph[neighbor].h = 50*rh_map.cost[index];
                            Graph[neighbor].parent = curr;
                            open_list.push({Graph[neighbor].g+Graph[neighbor].h, neighbor});
                        }
                    }
                    // cout<<"neighbor "<<neighbor<<endl<<"parent "<<curr<<endl;
                    // cout<<Heuristic_Map[{neighbor[0], neighbor[1]}].count<<" "<<max_time-neighbor[2]<<endl;
                }
            }
        }
        
//        cout<<(terminal.find(curr)==terminal.end() && !open_list.empty())<<endl;
        if(isTerminal(curr, terminal, time_lookup, thresh) && !open_list.empty()){
            cout<<"Termination condition reached"<<endl;
        }
    }
    points.close();

    cout<<"Terminated"<<endl;
    cout<<"Reached "<<curr<<endl;
    auto temp = curr;
    vector<int> goal = {robotposeX, robotposeY, 0};
    cout<<"Goal is "<<goal<<endl;
    // backpropagate to find the complete trajectory
    while(curr!=goal){
        // cout<<curr<<endl;
        trajectory.push_back({curr[0],curr[1]});
        curr = Graph[curr].parent;
        // __libcpp_thread_sleep_for(std::chrono::milliseconds(1000));
    }
    trajectory.push_back({curr[0], curr[1]});

    // reverse the vector since we had backpropagated
    reverse(trajectory.begin(), trajectory.end());
    int diff = target_steps - trajectory.size();

    // concatenate states to the end of trajectory if robot reaches earlier than the target
    if(diff>0){
        cout<<"Concatenating to planned trajecotry"<<endl;
        diff+=10;
        while(diff-- >= 0)
            trajectory.push_back({temp[0], temp[1]});
    }

    now = chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
    cout<<"Time taken to complete planning: "<<duration.count()<<endl;
    // points.close();
    created = true;
    action_ptr[0] = trajectory[current].first;
    action_ptr[1] = trajectory[current].second;
    current++;
    // cout.close();
    return;
}


// prhs contains input parameters (4):
// 1st is matrix with all the obstacles
// 2nd is a row vector <x,y> for the robot position
// 3rd is a matrix with the target trajectory
// 4th is an integer C, the collision threshold for the map
// plhs should contain output parameters (1):
// 1st is a row vector <dx,dy> which corresponds to the action that the robot should make
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    
    /* Check for proper number of arguments */
    if (nrhs != 6) {
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Six input arguments required.");
    } else if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required.");
    }
    
    /* get the dimensions of the map and the map matrix itself*/
    int x_size = mxGetM(MAP_IN);
    int y_size = mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the dimensions of the robotpose and the robotpose itself*/
    int robotpose_M = mxGetM(ROBOT_IN);
    int robotpose_N = mxGetN(ROBOT_IN);
    if(robotpose_M != 1 || robotpose_N != 2){
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidrobotpose",
                "robotpose vector should be 1 by 2.");
    }
    double* robotposeV = mxGetPr(ROBOT_IN);
    int robotposeX = (int)robotposeV[0];
    int robotposeY = (int)robotposeV[1];
    
    /* get the dimensions of the goalpose and the goalpose itself*/
    int targettraj_M = mxGetM(TARGET_TRAJ);
    int targettraj_N = mxGetN(TARGET_TRAJ);
    
    if(targettraj_M < 1 || targettraj_N != 2)
    {
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidtargettraj",
                "targettraj vector should be M by 2.");
    }
    double* targettrajV = mxGetPr(TARGET_TRAJ);
    int target_steps = targettraj_M;
    
    /* get the current position of the target*/
    int targetpose_M = mxGetM(TARGET_POS);
    int targetpose_N = mxGetN(TARGET_POS);
    if(targetpose_M != 1 || targetpose_N != 2){
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidtargetpose",
                "targetpose vector should be 1 by 2.");
    }
    double* targetposeV = mxGetPr(TARGET_POS);
    int targetposeX = (int)targetposeV[0];
    int targetposeY = (int)targetposeV[1];
    
    /* get the current timestep the target is at*/
    int curr_time = mxGetScalar(CURR_TIME);
    
    /* Create a matrix for the return action */ 
    ACTION_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)2, mxDOUBLE_CLASS, mxREAL); 
    double* action_ptr = (double*) mxGetData(ACTION_OUT);
    
    /* Get collision threshold for problem */
    int collision_thresh = (int) mxGetScalar(COLLISION_THRESH);
    
    /* Do the actual planning in a subroutine */
    planner(map, collision_thresh, x_size, y_size, robotposeX, robotposeY, target_steps, targettrajV, targetposeX, targetposeY, curr_time, &action_ptr[0]);
    // printf("DONE PLANNING!\n");
    return;   
}
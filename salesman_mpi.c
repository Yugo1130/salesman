// mpicc salesman_mpi.c -o salesman_mpi
// time mpiexec -np 4 ./salesman_mpi TSPmaps/gr24.tsp

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>

int n;  // 都市数
int **w;  // データ格納用

int ub = INT_MAX;  // 上界（これより悪くならない）
int lb = 0;  // 下界（これより良くならない）
int *best_route;  // 最短パスを記録
int *current_route;  // 現在地点までのパスを記録
bool *visited;  // 訪れたか

enum {
    TAG_REQ  = 10,   // タスク要求
    TAG_TASK = 20,   // タスク配布
    TAG_UB   = 30,   // UB更新
    TAG_TERM = 40    // 終了通知
};

// タスク管理用
typedef struct {
    int step;
    int cost;
    int current_node;
    int *route;
} Task;

Task *tasks = NULL;
int task_size = 0;

int read_file(const char *file_path){
    FILE *fp = fopen(file_path, "r");

    if (!fp) {
        perror("fopen");
        return 1;
    }

    // 都市数読み込み
    fscanf(fp, "%d", &n);

    // メモリ確保
    w = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; ++i) {
        w[i] = (int*)malloc(n * sizeof(int));
    }

    // 行列に距離を書く
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            int x;
            fscanf(fp, "%d", &x);
            w[i][j] = x;
            w[j][i] = x;
        }
    }

    fclose(fp);
    return 0;
}

// 貪欲ー法（最近傍法）を用いて初期ルートを導出
void greedy(int *route){
    bool *visited_greedy = calloc(n, sizeof(bool));  // 訪れたか
    int current_node = 0;  // 現在地点

    visited_greedy[0] = true; 
    route[0] = 0;
    
    for(int step = 1; step < n; step++){
        int best_cost = INT_MAX;
        int next_node = -1;
        // 現在地から一番近い都市を探す
        for(int candidate = 0; candidate < n; candidate++){
            if(!visited_greedy[candidate] && w[current_node][candidate] < best_cost){
                best_cost = w[current_node][candidate];
                next_node = candidate;
            }
        }
        // root[step] = next_node;
        visited_greedy[next_node] = true;
        route[step] = next_node;
        current_node = next_node;
    }
    free(visited_greedy);
}

int calc_cost(const int *route){
    int cost = 0;
    for(int i = 0; i < n; i++){
        cost += w[route[i]][route[(i + 1) % n]];
    }
    return cost;
}

void swap_edge(int *route, int l, int r){
    while(l < r){
        int temp = route[l];
        route[l] = route[r];
        route[r] = temp;
        l++;
        r--;
    }
}

// 2-opt法
int two_opt(int *route){
    int best_cost = calc_cost(route);
    bool improved = true;

    while(improved){
        improved = false;
        int best_diff = 0;
        int best_i;
        int best_j;
        for(int i = 0; i < n; i++){
            int a = route[i];
            int b = route[(i + 1) % n];
            for(int j = i + 2; j < n; j++){
                int c = route[j];
                int d = route[(j + 1) % n];

                int diff =  (w[a][b] + w[c][d]) - (w[a][c] + w[b][d]);  // 正なら改善
                if(best_diff < diff){
                    best_diff = diff;
                    best_i = i;
                    best_j = j;
                    improved = true;
                }
            }
        }

        if(improved){
            swap_edge(route, best_i + 1, best_j);
            best_cost -= best_diff;
        }
    }
    return best_cost;
}

// UBを計算
int calc_ub(){
    // 貪欲法（最近傍法）で初期解を導出し
    // さらに2-opt法を適用することで，近似解を導出する．
    greedy(best_route);
    int cost = two_opt(best_route);

    return cost;
}

// 訪問済の端点である0とcurrent_nodeから未訪問集合につなぐ際の最小コストを導出
int conect_u(int current_node, const bool *visited){
    int best_0 = INT_MAX;  // 端点0から未訪問集合U内の最短コスト
    int best_current = INT_MAX;  // 端点currentから未訪問集合U内の最短コスト

    for(int i = 1; i < n; i++){  // 端点0は除く
        if(visited[i]) continue;

        // 端点（0とcurrent_node，iは訪問先）
        if(w[0][i] < best_0) {
            best_0 = w[0][i];
        }
        if(w[current_node][i] < best_current) {
            best_current = w[current_node][i];
        }
    }

    return best_0 + best_current;
}

int min_two(int current_node, const bool *visited){
    int unvisited_cost = 0;  // コスト計算用
    int best_a;  // 常にbest_a <= best_b
    int best_b;
    int distance;

    for(int i = 1; i < n; i++){  // 端点0は除く
        if(visited[i]) continue;

        // 未訪問ノード（iは訪問元）
        best_a = INT_MAX;
        best_b = INT_MAX;
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            if(visited[j] && j != 0 && j != current_node) continue; // 0, current_nodeを除いた訪れているノードの場合

            distance = w[i][j];

            if(distance < best_a){
                best_b = best_a;
                best_a = distance;
            }
            else if(distance < best_b){
                best_b = distance;
            }
        }
        unvisited_cost += best_a + best_b;
    }
    return unvisited_cost;
}

int mst(const bool *visited){
    int cost = 0;

    // 未訪問集合Uの数を導出
    int m = 0;
    for(int i = 1; i < n; i++){
        if(!visited[i]) m++;
    }

    // 未訪問集合Uを作成
    int *U = malloc((m) * sizeof(int));
    int p = 0;
    for(int i = 1; i < n; i++){
        if(!visited[i]){
            U[p++] = i;
        }
    }

    // MSTにu[v]をつなげる最安コストを保持
    int *key = malloc((m) * sizeof(int));
    for (int i = 0; i < m; i++){
        key[i] = INT_MAX;
    }
    key[0] = 0; // u[0]から開始するため，u[0]のコストを0に設定

    bool *visited_mst = calloc(m, sizeof(bool));

    for(int i = 0; i < m; i++){
        // まだ訪れていないノードの中でkeyが最小のノードuを選択
        int u = -1;
        int best_node = INT_MAX;
        for(int j = 0; j < m; j++){
            if(visited_mst[j]) continue;
            if (key[j] < best_node) {
                best_node = key[j]; 
                u = j;
            }
        }
        visited_mst[u] = true;
        cost += key[u];

        // ノードuが追加されたため，各ノードのkeyを再計算
        for(int j = 0; j < m; j++){
            if(visited_mst[j]) continue;
            int diff = w[U[u]][U[j]];
            if(diff < key[j]){
                key[j] = diff;
            }
        }
    }
    free(U);
    free(visited_mst);
    free(key);
    return cost;
}

// LBを計算
int calc_lb(int cost, int current_node, const bool *visited){
    int hasU = false;
    for(int i = 1; i < n; i++){
        if(!visited[i]){
            hasU = true;
            break;
        }
    }
    // すべて訪問済みならノード0に戻る
    if(!hasU){
        return cost + w[current_node][0];
    }

    // 訪問済ノードの端点である0とcurrentから未訪問集合Uへの最短距離1本のコスト
    int cost_conect_u = conect_u(current_node, visited);

    // 未訪問ノードからUもしくは0もしくはcurrentへの最短距離2本のコストの平均を取ったものの合計
    int cost_min_two = (min_two(current_node, visited) + cost_conect_u) / 2;
    // 最小全域木を作成し，コストを導出
    int cost_mst = mst(visited) + cost_conect_u;

    // 各未訪問ノードに対して最短の2本を取る下界とMSTで得られた下界の最大値をLBとして使用
    if(cost_min_two > cost_mst){
        cost += cost_min_two;
    }
    else{
        cost += cost_mst;
    }

    return cost;
}

// 深さ優先探索
void dfs(int step, int cost, int current_node) {
    MPI_Status status;
    MPI_Request request;

    // rank0から最新のubがsendされていた場合は更新
    int latest_ub = INT_MAX;
    int flag;
    MPI_Iprobe(0, TAG_UB, MPI_COMM_WORLD, &flag, &status);
    if(flag){
        MPI_Recv(&latest_ub, 1, MPI_INT, 0, TAG_UB, MPI_COMM_WORLD, &status);
        if(latest_ub < ub){
            ub = latest_ub;
        }
    }
    
    // 枝切り
    if (cost >= ub) {
        return;
    }
    
    // lbを計算
    lb = calc_lb(cost, current_node, visited);
    
    
    // 枝切り
    if (lb >= ub) {
        return;
    }

    // 葉まで到達したら
    if (step == n) {
        // ノード0に戻る
        int total_cost = cost + w[current_node][0];
        if (total_cost < ub){
            ub = total_cost;
            printf("updatecost: %d\n", ub);
            int *buf = (int*)malloc((n + 1) * sizeof(int));
            buf[0] = ub;
            memcpy(buf + 1, current_route, n * sizeof(int));
            // ubが更新されたらrank0に送信
            MPI_Send(buf, n + 1, MPI_INT, 0, TAG_UB, MPI_COMM_WORLD);
            free(buf);
        }
        return;
    }

    // 全ルート探索
    for (int next_node = 1; next_node < n; next_node++) {
        if(!visited[next_node]){
            visited[next_node] = true;
            current_route[step] = next_node;
            dfs(step + 1, cost + w[current_node][next_node], next_node);
            visited[next_node] = false;
        }
    }
}

void push_task(int step, int cost, int current_node, int *route){
    tasks[task_size].step = step;
    tasks[task_size].cost = cost;
    tasks[task_size].current_node = current_node;
    tasks[task_size].route = (int*)malloc((step) * sizeof(int));
    memcpy(tasks[task_size].route, route, step * sizeof(int));
    task_size++;
}

void seed_dfs(int step, int cost, int current_node, int split_depth){
    // 枝切り
    if (cost >= ub) {
        return;
    }

    // lbを計算
    lb = calc_lb(cost, current_node, visited);

    // 枝切り
    if (lb >= ub) {
        return;
    }

    // 事前に設定した深さまで探索
    if (step == split_depth) {
        push_task(step, cost, current_node, current_route);
        return;
    }

    // 全ルート探索
    for (int next_node = 1; next_node < n; next_node++) {
        if(!visited[next_node]){
            visited[next_node] = true;
            current_route[step] = next_node;
            seed_dfs(step + 1, cost + w[current_node][next_node], next_node, split_depth);
            visited[next_node] = false;
        }
    }
}

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        fprintf(stderr, "引数にファイルパスを渡してください．");
        return 1;
    }

    if (size < 2){
        fprintf(stderr, "プロセス数を2以上にしてください．");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    read_file(argv[1]);

    best_route = (int*)malloc((n) * sizeof(int));
    current_route = (int*)malloc((n) * sizeof(int));
    visited = calloc(n, sizeof(bool));  // 訪れたかどうか

    ub = calc_ub();  // 上界を貪欲法によって計算

    if(rank == 0){
        printf("start ub: %d\n", ub);
    }

    int num_tasks = 1;
    int split_depth = 1; 

    while(num_tasks < (size - 1) * 20 && split_depth < n - 1){
        num_tasks *= (n - split_depth);
        split_depth++;
    }

    visited[0] = true;
    current_route[0] = 0;
    tasks = (Task*)malloc((num_tasks) * sizeof(Task));
    seed_dfs(1, 0, 0, split_depth);

    if(rank == 0){
        int task_num = 0;
        int flag;
        int num_term_processes = 0;
        bool processing = true;
        int temp;
        
        while(processing){
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if(flag){
                switch(status.MPI_TAG){
                    case TAG_REQ:{  // タスク要求
                        MPI_Recv(&temp, 1, MPI_INT, status.MPI_SOURCE, TAG_REQ, MPI_COMM_WORLD, &status);
                        if(task_num < task_size){  // タスク（番号）を分配
                            MPI_Send(&task_num, 1, MPI_INT, status.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD);
                            task_num++;
                        } else{  // 処理すべきタスクが残っていない場合は終了通知を送信
                            MPI_Send(&temp, 1, MPI_INT, status.MPI_SOURCE, TAG_TERM, MPI_COMM_WORLD);
                            num_term_processes++;
                        }
                        if(num_term_processes >= size - 1){  // すべてのプロセスに終了通知を送ったらwhileから抜ける
                            processing = false;
                        }
                        // printf("rank0: %d\n", num_term_processes);
                        break;
                    }
                    case TAG_UB:{ // UB更新
                        int *buf = (int*)malloc((n+1) * sizeof(int));
                        MPI_Recv(buf, n + 1, MPI_INT, MPI_ANY_SOURCE, TAG_UB, MPI_COMM_WORLD, &status);
                        int latest_ub = buf[0];
                        if(latest_ub < ub){
                            ub = latest_ub;
                            // すべてのrankにお知らせ
                            for(int dest = 1; dest < size; dest++){
                                MPI_Isend(&ub, 1, MPI_INT, dest, TAG_UB, MPI_COMM_WORLD, &request);
                            }
                            memcpy(best_route, &buf[1], n * sizeof(int));
                        }
                        break;
                    }
                }
            }

        }
        
    } else {
        int flag;
        bool processing = true;
        int task_num;
        int temp;
        // 初回のタスク要求
        MPI_Send(&temp, 1, MPI_INT, 0, TAG_REQ, MPI_COMM_WORLD);

        while(processing){
            MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if(flag){
                switch(status.MPI_TAG){
                    case TAG_TASK:{  // タスク受信
                        MPI_Recv(&task_num, 1, MPI_INT, 0, TAG_TASK, MPI_COMM_WORLD, &status);
                        memset(visited, 0, n * sizeof(bool));
                        for (int i = 0; i < tasks[task_num].step; ++i) {
                            current_route[i] = tasks[task_num].route[i];
                            visited[current_route[i]] = true;
                        }
                        dfs(tasks[task_num].step, tasks[task_num].cost, tasks[task_num].current_node);
                        // 処理が終わったら新たなタスクを要求
                        MPI_Send(&temp, 1, MPI_INT, 0, TAG_REQ, MPI_COMM_WORLD);
                        break;
                    }
                    case TAG_TERM:{
                        MPI_Recv(&temp, 1, MPI_INT, 0, TAG_TERM, MPI_COMM_WORLD, &status);
                        processing = false;
                        break;
                    }
                }
            }
        }
    }

    if(rank == 0){
        printf("bestcost: %d\n", ub);

        for (int i = 0; i < n; i++) {
            printf("%d ", best_route[i]);
        }
        printf("0\n");
    }

    for (int i = 0; i < n; ++i) free(w[i]);
    free(w);
    free(best_route);
    free(current_route);
    for (int i = 0; i < task_size; ++i) free(tasks[i].route);
    free(tasks);
    MPI_Finalize();
    return 0;
}

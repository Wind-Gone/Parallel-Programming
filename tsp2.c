#include <mpi.h>
#include <bits/stdc++.h>
using namespace std;
#define N_OF_CITY 15           //定义城市数量
int cost[N_OF_CITY][N_OF_CITY]; //定义城市之间的路径代价

#define TAG 1         //用于通信的标识
#define RESULT 2      //用于通信的标识
#define Available -1  //可用
#define Unavailable 0 //路径不可用

//定义结构体
typedef struct
{
    int path[N_OF_CITY + 2]; //完整路径
    int city_count;          //已用城市数量
    int curcost;             //当前路径代价
    string city_name;        //城市名（未用到）
} Travel;

// 初始化路径
void initTravel(Travel *tour)
{
    for (int i = 0; i < N_OF_CITY + 2; i++)
    {
        tour->path[i] = Available;
    }
    tour->path[0] = 0;
    tour->city_count = 1;
    tour->curcost = 0;
}

// 增加一座城市
void pushTravel(Travel *tour, int nbr)
{
    int ini = tour->path[(tour->city_count) - 1];
    tour->curcost += cost[ini][nbr];
    tour->path[tour->city_count] = nbr;
    tour->city_count++;
}

// 移出一座城市
void removeTravel(Travel *tour)
{
    int src = tour->path[(tour->city_count) - 1];
    int dst = tour->path[(tour->city_count) - 2];
    tour->curcost -= cost[dst][src];
    tour->path[(tour->city_count) - 1] = -1;
    tour->city_count--;
}

//城市路径覆盖
void copyTravel(Travel *copy, Travel *src)
{
    for (int i = 0; i < N_OF_CITY + 2; i++)
    {
        copy->path[i] = src->path[i];
    }
    copy->city_count = src->city_count;
    copy->curcost = src->curcost;
}

//是否是最佳路径
int Best_tour(Travel *tour, Travel *bestTravel)
{
    pushTravel(tour, 0);
    return tour->curcost < bestTravel->curcost;
}

//城市是否可行（即可能获得更小的回路代价）
int Feasible(Travel curr_tour, int nbr, Travel bestTravel)
{
    if (curr_tour.curcost > bestTravel.curcost)
        return 0;
    for (int i = 0; i < curr_tour.city_count; i++)
    {
        if (curr_tour.path[i] == nbr)
            return 0;
    }
    return 1;
}

//打印最终结果
void print_cost(int path_res[], int cost)
{
    cout << "代价：" << cost << endl
         << "路径为：";
    for (int i = 0; i <= N_OF_CITY; i++)
    {
        cout << path_res[i] << " ";
    }
    cout << endl;
}

//计算最终的路径代价
int cal_cost(int path[])
{
    int rst = 0;
    for (int i = 0; i < N_OF_CITY; i++)
    {
        rst += cost[path[i]][path[i + 1]];
    }
    return rst;
}

//更新路径
void update_path(int src[], int dst[])
{
    for (int i = 0; i < N_OF_CITY + 2; i++)
        dst[i] = src[i];
}

//初始化代价
void fill_distance()
{
    //每座城市之间代价初始化
    srand(1);
    for (int i = 0; i < N_OF_CITY; i++)
    {
        for (int j = 0; j < N_OF_CITY; j++)
        {
            cost[i][j] = rand() % 8 + 1;
        }
    }
}

int main(int argc, char *argv[])
{
    int my_rank, comm_sz;
    fill_distance();
    double t1, t2 = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Status status;

    int path_cur[N_OF_CITY + 2];
    memset(path_cur, 0, sizeof(int));

    int send_buf_sz = (N_OF_CITY - 1) / comm_sz; //每个进程将要分到的分支
    int recv_buf_sz = send_buf_sz;               //发送的和收到的进程处理的分支个数应该相同
    Travel final_path;
    initTravel(&final_path);
    final_path.curcost = INT_MAX;

    int msg_avail = 0;
    int *send_buf_p;

    if (my_rank == 0)
    {
        t1 = MPI_Wtime(); //开始计时
        send_buf_p = (int *)malloc(comm_sz * send_buf_sz * sizeof(int));
        for (int i = 0; i < N_OF_CITY - 1; i++)
        {
            send_buf_p[i] = i + 1;
        }
        cout << "Computing solution using " << comm_sz << " processes" << endl;
        cout << "start now!" << endl;
    }
    int recv_buf_p[recv_buf_sz];
    //scatter 阻塞式的，0号进程分配任务
    MPI_Scatter(send_buf_p, send_buf_sz, MPI_INT, recv_buf_p, recv_buf_sz, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < send_buf_sz; i++)
    {
        printf("[%d]进程: 开始执行任务\n", my_rank);
        stack<Travel> stack;
        Travel bestTravel;
        initTravel(&bestTravel);
        bestTravel.curcost = INT_MAX;
        bestTravel.city_count = N_OF_CITY;

        Travel tmptravel;
        initTravel(&tmptravel);
        pushTravel(&tmptravel, recv_buf_p[i]);
        stack.push(tmptravel);
        while (!stack.empty())
        {

            Travel cur_tour = stack.top();
            stack.pop();
            if (cur_tour.city_count == N_OF_CITY)
            {
                if (Best_tour(&cur_tour, &bestTravel))
                {
                    bestTravel = cur_tour;
                }
            }
            else
            {

                for (int nbr = N_OF_CITY - 1; nbr >= 1; nbr--)
                {
                    if (Feasible(cur_tour, nbr, bestTravel))
                    {
                        pushTravel(&cur_tour, nbr);
                        stack.push(cur_tour);
                        removeTravel(&cur_tour);
                    }
                }
            }
        }
        //将结果赋给最终路径（前提是有程序运行正常，有最佳路径）
        final_path = final_path.curcost > bestTravel.curcost ? bestTravel : final_path;
    }

    //等待所有进程结束

    MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d]进程: 终于结束了\n", my_rank);
    int final_cost;
    int final_best_path[100];
    final_cost = final_path.curcost;
    final_best_path[N_OF_CITY + 2] = final_path.curcost;
    MPI_Reduce(&final_path.curcost, &final_cost, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD); //全局规约求最小值

    if (my_rank == 0)
        print_cost(final_best_path, final_cost);
    MPI_Finalize();
    t2 = MPI_Wtime(); //计时结束
    if (my_rank == 0)
        cout << "计算时间为 " << t2 - t1 << "s" << endl;
    return 0;
}

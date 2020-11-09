#include <mpi.h>
#include <bits/stdc++.h>
using namespace std;
#define N_OF_CITY 15  //定义城市数量
int cost[N_OF_CITY][N_OF_CITY]; //定义城市之间的路径代价

#define TAG 1  //用于通信的标识
#define Available -1 //可用
#define Unavailable 0 //路径不可用

//定义结构体
typedef struct
{
    int path[N_OF_CITY + 2]; //完整路径
    int city_count; //已用城市数量
    int curcost; //当前路径代价
    string city_name;  //城市名（未用到）
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
void print_cost(Travel tour)
{
    cout << endl;
    cout<<"代价："<<tour.curcost<<endl<<"路径为：";
    for (int i = 0; i < tour.city_count; i++)
    {
        cout << tour.path[i] << " ";
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

    int send_buf_sz = (N_OF_CITY - 1) / comm_sz;//每个进程将要分到的分支
    int recv_buf_sz = send_buf_sz; //发送的和收到的进程处理的分支个数应该相同
    Travel final_path;
    initTravel(&final_path);
    final_path.curcost = INT_MAX;

    int msg_avail = 0;
    int *send_buf_p;

    if (my_rank == 0)
    {
        t1 = MPI_Wtime();//开始计时
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
            //如果运行过程中发现某个线程发来最佳路径，该线程就接收并更新自己的最优路径
            MPI_Iprobe(MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &msg_avail, &status);
            while (msg_avail)
            {
                MPI_Recv(&path_cur, N_OF_CITY + 2, MPI_INT, status.MPI_SOURCE, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (cal_cost(path_cur) < bestTravel.curcost)
                {
                    update_path(path_cur, bestTravel.path);
                    bestTravel.curcost = cal_cost(path_cur);
                }
                MPI_Iprobe(MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &msg_avail, &status);
            }

            Travel cur_tour = stack.top();
            stack.pop();
            if (cur_tour.city_count == N_OF_CITY)
            //如果所有城市都走完了，就需要更新自己的最佳路径
            {
                if (Best_tour(&cur_tour, &bestTravel))
                {
                    bestTravel = cur_tour;
                    for (int i = 0; i < comm_sz; i++)
                    {
                        if (i != my_rank)
                        {
                            update_path(bestTravel.path, path_cur);
                            MPI_Send(&path_cur, N_OF_CITY + 2, MPI_INT, i, TAG, MPI_COMM_WORLD);
                        }
                    }
                }
            }
            else
            //如果城市没走完，就前往下一个城市
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
        //将结果赋给每个线程的最终路径（前提是有程序运行正常，有最佳路径）
        final_path = final_path.curcost > bestTravel.curcost ? bestTravel : final_path;
    }

    // 这里为了防止0号进程没有运行速度过快没有收到最佳路径，我们再执行一次发送接受各自最佳路径和全局规约操作
    if(my_rank!=0){
		update_path(final_path.path,path_cur);
		MPI_Send(&path_cur, N_OF_CITY+2, MPI_INT, 0, TAG, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d]进程: 终于结束了\n", my_rank);
	if(my_rank==0){
		MPI_Iprobe( MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,&msg_avail,&status);
		while(msg_avail){
			MPI_Recv(&path_cur, N_OF_CITY+2, MPI_INT, status.MPI_SOURCE,TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(cal_cost(path_cur)<final_path.curcost){
				update_path(path_cur,final_path.path);
				final_path.curcost=cal_cost(path_cur);
			}	
			MPI_Iprobe( MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,&msg_avail,&status);
		}
         //打印最佳路径和代价
        print_cost(final_path );
	}
        
    MPI_Finalize();
    t2 = MPI_Wtime(); //计时结束
    if (my_rank == 0)
        cout << "计算时间为 " << double(t2 - t1) << "s" << endl;
    return 0;
}

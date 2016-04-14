/*构建fat-tree和流量矩阵*/

#include<iostream>
#include<Eigen/Dense>
#include<fstream>
#include<vector>
#include<random>
#include<chrono>
using namespace std;
using namespace Eigen;

/* Core层有(k^2/4)个交换机
*  Aggregation 层有(k^2/2)个交换机
*  Edge层有(k^2/2)个交换机
*  host共有(k^3/4)个（暂时不靠考虑host）
*  每个pod里面有k/2个交换机
*/

/* k = 4时，交换机共有20个
*  k = 8时，交换机共有80个
*/

/*
* 交换机之间如何连接？(从上往下，从左往右编号)
* 1: core和aggregation：每个aggregation switch（即为as）有k个port（其实每个交换机都这样）,
*    as的up port有k/2个，则每个as和k/2个cs(core switch)连接。从as的角度出发，每个pod中的as按顺序选k/2个cs
*    比如，第0个as连接编号为0~k/2-1的cs，第1个as连接k/2~2*(k/2)-1的cs，第2个as连接2*(k/2)~3*(k/2)-1的cs
*    最后，第(k/2-1)个as连接编号为(k/2-1)*(k/2)~k/2*(k/2)-1的cs
* 2：pod i对应的as为：k*k/4+k/2*i ~ k*k/4+k/2*(i+1)-1
* 2：pod i对应的es(Edge switch)为：3*k*k/4+k/2*i ~ 3*k*k/4+k/2*(i+1)-1
*/

const int inf = 1e9;
const int MAX_SWITCH = 100;
const int MAX_ROUTING_CLASS = 100;
class Prob
{
private:
    int k;// k pods
    int num_core_switch;
    int num_aggregation_switch;
    int num_edge_switch;
    int num_link;
    int num_flow;

    int num_routing_class;
    //int num_host;
    MatrixXd topo;//拓扑矩阵
    MatrixXd capacity;
    MatrixXd link_flow;
    MatrixXd switch_to_link;//两台交换机之间的链路编号
    /*
    * 设置traffic demand matrices
    * 暂时只考虑Edge层之间的交换机互相发送flows
    * 假设每个值是均匀分布于[0.01,1]
    */
    MatrixXd traffic_demand_matrices[MAX_ROUTING_CLASS];
    MatrixXd TM;//TM 为 routing class 列的矩阵,行数为: num_edge_switch*num_edge_switch
    int N;//交换机总数
    unsigned seed;
public:
    friend class Run;
    Prob(int K);
    void generate_topo();
    void generate_VLAN();
    void generate_link_flow();
    void generate_TD();//生成 traffic demand matrices
    MatrixXd get_topo();
    MatrixXd get_capacity();
    //Dojkstra: 找到从s出发到所有交换机的最短路径的前驱节点
    void Dijkstra(vector<vector<int>> pre[],int s,MatrixXd G);
    //get_path：随机生成一个路径
    void get_path(vector<int> &path,vector<vector<int>> pre[],int s,int d);
};

Prob::Prob(int K=4)
{
    seed = chrono::system_clock::now().time_since_epoch().count();

    k = K;
    num_core_switch = k*k/4;
    num_aggregation_switch = k*k/2;
    num_edge_switch = k*k/2;
    num_link = k*k*k;
    num_flow = (k*k/2)*(k*k/2);//只考虑edge switch交换机互相之间发送flow
    num_routing_class = k;//默认k个路由类
    N = num_core_switch+num_aggregation_switch+num_edge_switch;//交换机总数
    //生成拓扑矩阵和链路带宽
    generate_topo();
    generate_link_flow();
    generate_TD();

}
void Prob::generate_topo()
{


    topo = MatrixXd::Zero(N,N);
    capacity = MatrixXd::Zero(num_link,1);
    switch_to_link = MatrixXd::Zero(N,N);
    int id=0;
    //汇聚层和核心层连接
    //每个pod:i对应的汇聚层交换机编号为：k*k/4+k/2*i ~ k*k/4+k/2*(i+1)-1
    for(int i=0; i<k; i++)
    {
        for(int sw = k*k/4+k/2*i; sw<k*k/4+k/2*(i+1); sw++)
        {
            int sw_in_pod = sw%(k/2);//0~k/2-1
            // sw_in_pod 和 k/2个交换机连接
            for(int z=k/2*sw_in_pod; z<k/2*(sw_in_pod+1); z++) //z为核心层交换机
            {
                topo(sw,z) = 1;
                topo(z,sw) = 1;
                switch_to_link(sw,z) = id;
                switch_to_link(z,sw) = id+1;
                capacity(id,0) = 10000;
                capacity(id+1,0) = 10000;
                id += 2;

            }
        }
    }
    //汇聚层和Edge层连接
    for(int i=0; i<k; i++)
    {
        for(int as = k*k/4+k/2*i; as<k*k/4+(k/2)*(i+1); as++)
        {
            for(int es = 3*k*k/4+k/2*i; es<3*k*k/4+(k/2)*(i+1); es++)
            {
                topo(as,es)=1;
                topo(es,as)=1;
//                ofstream out0;
//                out0.open("F:\\path-index.txt",ios::app);
//                out0<<as<<"-"<<es<<": "<<id<<endl;
//                out0<<es<<"-"<<as<<": "<<id+1<<endl;
//                out0.close();
                switch_to_link(as,es) = id;
                switch_to_link(es,as) = id+1;
                capacity(id,0) = 10000;
                capacity(id+1,0) = 10000;
                id += 2;
            }
        }
    }

}
MatrixXd Prob::get_topo()
{
    return topo;
}
MatrixXd Prob::get_capacity()
{
    return capacity;
}


void Prob::generate_link_flow()
{
    /*ECMP 不需要考虑VLAN，直接从多个最短路径（equal cost）中随机选择一个*/
    /* 先只考虑一种routing class，即es之间只互相发送一遍 */
    link_flow = MatrixXd::Zero(num_link,num_edge_switch*num_edge_switch*num_routing_class);
    //ofstream out;
    //out.open("F:\\test_path.txt");
    vector< vector<int> > pre[MAX_SWITCH];

    for(int i=0; i<N; i++)
    {
        vector<int> temp;
        temp.clear();
        for(int j=0; j<N; j++)
        {
            pre[i].push_back(temp);
        }
        pre[i][i].push_back(i);
        //cout<<"hi"<<endl;
    }
    for(int s=0; s<N; s++)
    {
        Dijkstra(pre,s,topo);
    }

//    for(int i=0;i<N;i++){
//        for(int j=0;j<N;j++){
//            if(i!=j){
//                if(pre[i][j].size() == 1 && pre[i][j][0] == i)
//                    pre[i][j][0] = j;
//            }
//
//        }
//    }

//    for(int i=0;i<N;i++){
//
//        for(int j=0;j<N;j++){
//            out<<i<<" -> "<<j<<" : "<<endl;
//            for(int z=0;z<pre[i][j].size();z++)
//                out<<pre[i][j][z]<<" ";
//            out<<endl;
//        }
//
//    }

//    vector<int> path;
//    get_path(path,pre,17,19);
//    reverse(path.begin(),path.end());
//    cout<<"17 -> 19"<<endl;
//    for(int i=0;i<path.size();i++)
//        cout<<path[i]<<" ";
//    cout<<endl;
    for(int i=0; i<num_routing_class; i++)
    {
        int len = num_edge_switch*num_edge_switch;
        MatrixXd temp_link_flow;
        temp_link_flow = MatrixXd::Zero(num_link,len);

        for(int s = 3*k*k/4; s<=5*k*k/4-1; s++)
        {
            for(int d = 3*k*k/4; d<=5*k*k/4-1; d++)
            {
                if(s == d)
                    continue;
                vector<int> path;
                path.clear();
                get_path(path,pre,s,d);
                reverse(path.begin(),path.end());
//                cout<<"17 -> 19"<<endl;
//                for(int i=0; i<path.size(); i++)
//                    cout<<path[i]<<" ";
//                cout<<endl;

                //先把交换机编号转化为从0开始的序号，要不然矩阵下标会溢出
                int source = s - 3*k*k/4;
                int destination = d - 3*k*k/4;

                int column = source*num_edge_switch + destination;
               // cout<<"i = "<<i<<" column = "<<column<<endl;
                int len = path.size();
                for(int t=0; t<len/2; t++)
                {
                    //cout<<"t = "<<t<<" len = "<<len<<endl;
                    int link = switch_to_link(path[2*t],path[2*t+1]);
                    //cout<<"t = "<<t<<" "<<link<<" "<<column<<endl;
                    temp_link_flow(link,column) = 1;
                }

            }
        }
        for(int row=0; row<num_link; row++)
        {
            for(int col=i*num_flow; col<(i+1)*num_flow; col++)
            {
                link_flow(row,col) = temp_link_flow(row,col-i*num_flow);
            }
        }

    }
}

void Prob::generate_TD()
{
    //unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> dis(0.01,1);
    for(int i=0; i<num_routing_class; i++)
    {
        traffic_demand_matrices[i] = MatrixXd::Zero(num_edge_switch,num_edge_switch);
        for(int row = 0; row<num_edge_switch; row++)
        {
            for(int col=0; col<num_edge_switch; col++)
            {
                if(row != col)
                {
                    traffic_demand_matrices[i](row,col) = dis(generator);
                }
            }
        }
    }
    TM = MatrixXd::Zero(num_edge_switch*num_edge_switch,num_routing_class);
    for(int c=0; c<num_routing_class; c++)
    {
        MatrixXd tm = traffic_demand_matrices[c];
        tm.resize(num_edge_switch*num_edge_switch,1);
        for(int r = 0; r<num_edge_switch*num_edge_switch; r++)
        {
            TM(r,c) = tm(r,0);
        }
    }
}
// 生成从s出发的最短路径
void Prob::Dijkstra(vector<vector<int>> pre[],int s,MatrixXd G)
{

    int dis[MAX_SWITCH] = {0};
    for(int i=0; i<N; i++)//N个交换机
    {
        if(i == s)
            continue;
        else
            dis[i] = inf;
    }
    bool is_visited[MAX_SWITCH]= {false};
    for(int j=0; j<N; j++)
    {
        int MIN = inf;
        int u = -1;
        for(int z=0; z<N; z++)
        {
            if(is_visited[z] == false && dis[z] < MIN)
            {
                MIN = dis[z];
                u = z;
            }
        }
        if(u == -1)
            break;
        is_visited[u] = true;

        for(int x=0; x<N; x++)
        {
            if(G(u,x) > 0)
                if(dis[u] + G(u,x) < dis[x])
                {
                    dis[x] = dis[u] + G(u,x);
                    pre[s][x].clear();
                    pre[s][x].push_back(u);

                }
                else if(dis[u] + G(u,x) == dis[x])
                {
                    pre[s][x].push_back(u);
                }
        }

    }

}

void Prob::get_path(vector<int> &path,vector<vector<int> > pre[],int s,int d)
{
    int len = pre[s][d].size();
    path.push_back(d);
    if(len == 1 && pre[s][d][0] == s)
    {
        path.push_back(pre[s][d][0]);
        return;
    }
    default_random_engine generator(seed);
    uniform_int_distribution<int> dis(0,len-1);
    int next = pre[s][d][dis(generator)];
    path.push_back(next);
    get_path(path,pre,s,next);
}

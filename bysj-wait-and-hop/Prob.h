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
const int MAX_VLAN = 20;
const int inf = 1e9;
const int MAX_SWITCH = 100;
const int MAX_ROUTING_CLASS = 100;
class Prob{
private:
    int k;// k pods
    int num_core_switch;
    int num_aggregation_switch;
    int num_edge_switch;
    int num_link;
    int num_flow;
    int num_vlan;
    int num_routing_class;
    //int num_host;
    MatrixXd topo;//拓扑矩阵
    MatrixXd capacity;
    MatrixXd VLAN[MAX_VLAN];
    MatrixXd link_flow[MAX_VLAN];//每个VLAN对应一个
    MatrixXd switch_to_link;//两台交换机之间的链路编号
    /*
    * 设置traffic demand matrices
    * 暂时只考虑Edge层之间的交换机互相发送flows
    * 假设每个值是均匀分布于[0.001,0.01]
    */
    MatrixXd traffic_demand_matrices[MAX_ROUTING_CLASS];
    MatrixXd TM;//TM 为 routing class 列的矩阵,行数为: num_edge_switch*num_edge_switch
    int N;//交换机总数
    bool is_sw_in_vlan[MAX_SWITCH][MAX_VLAN];
public:
    friend class Run;
    Prob(int K);
    void generate_topo();
    void generate_VLAN();
    void generate_link_flow();
    void generate_TD();//生成 traffic demand matrices
    MatrixXd get_topo();
    MatrixXd get_capacity();
    MatrixXd Dijkstra(MatrixXd G);//某个VLAN下的最短路径
    vector<int> get_shortest_path(MatrixXd next,int s,int d);
};

Prob::Prob(int K=4){
    k = K;
    num_vlan = k;
    num_core_switch = k*k/4;
    num_aggregation_switch = k*k/2;
    num_edge_switch = k*k/2;
    num_link = k*k*k;
    num_flow = (k*k/2)*(k*k/2);//只考虑edge switch交换机互相之间发送flow
    num_routing_class = k;//默认k个路由类
    N = num_core_switch+num_aggregation_switch+num_edge_switch;//交换机总数
    //生成拓扑矩阵和链路带宽
    generate_topo();
    //生成VLANs
    generate_VLAN();
    generate_link_flow();
    generate_TD();

}
void Prob::generate_topo(){


    topo = MatrixXd::Zero(N,N);
    capacity = MatrixXd::Zero(num_link,1);
    switch_to_link = MatrixXd::Zero(N,N);
    int id=0;
    //汇聚层和核心层连接
    //每个pod:i对应的汇聚层交换机编号为：k*k/4+k/2*i ~ k*k/4+k/2*(i+1)-1
    for(int i=0;i<k;i++){
        for(int sw = k*k/4+k/2*i;sw<k*k/4+k/2*(i+1);sw++){
            int sw_in_pod = sw%(k/2);//0~k/2-1
            // sw_in_pod 和 k/2个交换机连接
            for(int z=k/2*sw_in_pod;z<k/2*(sw_in_pod+1);z++){//z为核心层交换机
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
    for(int i=0;i<k;i++){
        for(int as = k*k/4+k/2*i;as<k*k/4+(k/2)*(i+1);as++){
            for(int es = 3*k*k/4+k/2*i;es<3*k*k/4+(k/2)*(i+1);es++){
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
MatrixXd Prob::get_topo(){
    return topo;
}
MatrixXd Prob::get_capacity(){
    return capacity;
}

void Prob::generate_VLAN(){
    for(int i=0;i<k;i++)
        VLAN[i] = MatrixXd::Zero(N,N);
    //每个core switch作为root
    for(int i=0;i<k;i++){
        //i表示core switch id
        //每个pod遍历一次
        is_sw_in_vlan[i][i] = true;
        for(int pod = 0;pod<k;pod++){
            for(int as = k*k/4+k/2*pod;as<k*k/4+(k/2)*(pod+1);as++){
                if(topo(as,i)>0){
                    VLAN[i](as,i) = 1;
                    VLAN[i](i,as) = 1;
                    is_sw_in_vlan[as][i] = true;
                    for(int es = 3*k*k/4+k/2*pod;es<3*k*k/4+(k/2)*(pod+1);es++){
                        VLAN[i](as,es) = 1;
                        VLAN[i](es,as) = 1;
                        is_sw_in_vlan[es][i] = true;
                    }
                }
            }
        }
    }

}
void Prob::generate_link_flow(){
    for(int i=0;i<k;i++){
        MatrixXd next;
        next = Dijkstra(VLAN[i]);
        //vector<int> path = get_shortest_path(next,s,d);
        MatrixXd temp_link_flow;


        temp_link_flow = MatrixXd::Zero(num_link,num_flow);
        link_flow[i] = MatrixXd::Zero(num_link,num_flow);
        //cout<<"i = "<<i<<" num_link = "<<num_link<<"  num_flow = "<<num_flow<<endl;
        for(int s = 3*k*k/4;s<=5*k*k/4-1;s++){
            for(int d = 3*k*k/4;d<=5*k*k/4-1;d++){
                if(s == d)
                    continue;
                vector<int> path = get_shortest_path(next,s,d);
                //先把交换机编号转化为从0开始的序号，要不然矩阵下标会溢出
                int source = s - 3*k*k/4;
                int destination = d - 3*k*k/4;

                int column = source*num_edge_switch + destination;
                //cout<<"i = "<<i<<" column = "<<column<<endl;
                int len = path.size();
                for(int t=0;t<len/2;t++){
                    //cout<<"t = "<<t<<" len = "<<len<<endl;
                    int link = switch_to_link(path[2*t],path[2*t+1]);
                    //cout<<"t = "<<t<<" "<<link<<" "<<column<<endl;
                    temp_link_flow(link,column) = 1;
                }

            }
        }
        link_flow[i] = temp_link_flow;
    }

}

void Prob::generate_TD(){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> dis(0.01,1);
    for(int i=0;i<num_routing_class;i++){
        traffic_demand_matrices[i] = MatrixXd::Zero(num_edge_switch,num_edge_switch);
        for(int row = 0;row<num_edge_switch;row++){
            for(int col=0;col<num_edge_switch;col++){
                if(row != col){
                    traffic_demand_matrices[i](row,col) = dis(generator);
                }
            }
        }
    }
    TM = MatrixXd::Zero(num_edge_switch*num_edge_switch,num_routing_class);
    for(int c=0;c<num_routing_class;c++){
        MatrixXd tm = traffic_demand_matrices[c];
        tm.resize(num_edge_switch*num_edge_switch,1);
        for(int r = 0;r<num_edge_switch*num_edge_switch;r++){
            TM(r,c) = tm(r,0);
        }
    }
}
MatrixXd Prob::Dijkstra(MatrixXd G){
    MatrixXd path;
    path = MatrixXd::Zero(N,N);
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            if(i == j)
                continue;
            else if(G(i,j) == 0)
                path(i,j) = -1;
        }

    int dis[MAX_SWITCH][MAX_SWITCH] = {0};
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(i == j)
                dis[i][j] = 0;
            else
                dis[i][j] = inf;
        }
    }
    //每个顶点都计算一次最短路径
    for(int i=0;i<N;i++){
        bool is_visited[MAX_SWITCH]={false};
        for(int j=0;j<N;j++){
            int MIN = inf;
            int u = -1;
            for(int z=0;z<N;z++){
                if(is_visited[z] == false && dis[i][z] < MIN){
                    MIN = dis[i][z];
                    u = z;
                }
            }
            if(u == -1)
                continue;
            is_visited[u] = true;
            //以上步骤找到离i最近的点u
            for(int x=0;x<N;x++){
                if(G(u,x) > 0 && dis[i][u] + G(u,x) < dis[i][x]){
                    dis[i][x] = dis[i][u] + G(u,x);
                    path(i,x) = u;
                    //cout<<"i = "<<i<<", x = "<<x<<",u = "<<u<<endl;
                }
            }
        }
    }
    return path;
}

vector<int> Prob::get_shortest_path(MatrixXd next,int s,int d){
    if(s == d){
        vector<int> t;
        t.clear();
        return t;
    }
    int k = next(s,d);
    if(k == s || k == d){
        vector<int> v;
        v.clear();
        v.push_back(s);
        v.push_back(d);
        return v;
    }else{
        vector<int> first = get_shortest_path(next,s,k);
        vector<int> second = get_shortest_path(next,k,d);
        vector<int> path;
        path.clear();
        for(int i=0;i<first.size();i++)
            path.push_back(first[i]);
        for(int i=0;i<second.size();i++)
            path.push_back(second[i]);
        return path;
    }
}

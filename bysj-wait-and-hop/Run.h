#include<vector>
#include<map>
#include<Eigen/Dense>
#include<algorithm>
#include<chrono>
#include<random>
#include<string>
using namespace Eigen;
using namespace std;

class Run{
private:
    double opt_len;
    int len;
    double total_time;
    double last_time;
    double beta;
    double alpha;
    vector<double> total_thr;
    vector<double> opt_time;
    Prob problem;
    vector<int> vlan_assignment;
    MatrixXd tm;
public:
    Run(Prob prob,double opt_length);
    vector<int> markov_vlan_assign(Prob problem);
    MatrixXd getRoute(Prob problem, vector<int> vlan_assignment);
    MatrixXd maxminRate(Prob problem,MatrixXd x,MatrixXd R);
    string vector2string(vector<int> x);
    double get_total_time();
    vector<double> get_total_thr();
    vector<double> get_opt_time();
};
Run::Run(Prob prob,double opt_length){
    problem = prob;
    opt_len = opt_length;
    len = problem.num_routing_class*problem.num_edge_switch*problem.num_edge_switch;
    total_time = 0;
    last_time = 0;
    beta = 0.1;
    alpha = 0.0001;
    tm = problem.TM;
    tm.resize(len,1);
    total_thr.clear();
    total_thr.push_back(0.0);
    opt_time.clear();
    opt_time.push_back(0.0);
    vlan_assignment.clear();
    vlan_assignment = markov_vlan_assign(problem);
    double max_link_u = -1;
//    for(int i=0;i<vlan_assignment.size();i++)
//        cout<<vlan_assignment[i]<<" ";
//    cout<<endl;
    while(1){
        if(opt_len<0){
            vlan_assignment = markov_vlan_assign(problem);

        }else if(total_time>opt_len + last_time){
            vlan_assignment = markov_vlan_assign(problem);
            last_time = total_time;
            cout<<"Hi"<<endl;
        }
        MatrixXd active = tm;
        for(int i=0;i<active.rows();i++){
            if(active(i,0)>0)
                active(i,0) = 1;
        }

        MatrixXd R;
        R = getRoute(problem,vlan_assignment);

        MatrixXd rate;
        rate = maxminRate(problem,active,R);
        //计算最大利用率
        MatrixXd links_utilization;
        links_utilization = R*rate;
        for(int i=0;i<problem.num_link;i++){
            links_utilization(i,0) = links_utilization(i,0)/problem.capacity(i,0);
            cout<<"链路利用率："<<links_utilization(i,0)<<endl;
        }

        //cout<<active.sum()<<" flows left\n";
        double duration = inf*1.0;
        bool changed = false;
        for(int i=0;i<len;i++){
            if(tm(i,0) > 0){
                double t = tm(i,0)/rate(i,0);
                if(0<t && t < duration){
                    duration = t;
                    changed = true;
                }
            }
        }
        if(changed == false)
            break;
        total_time += duration;
        for(int i=0;i<len;i++){
            if(tm(i,0) > 0)
                tm(i,0) = tm(i,0) - rate(i,0)*duration;
        }
        total_thr.push_back(rate.sum());
        opt_time.push_back(total_time);
    }
}
vector<int> Run::markov_vlan_assign(Prob problem){
    cout<<"vlan begin"<<endl;
    int numhash = problem.num_routing_class;
    int numvlan = problem.num_vlan;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    map<string,double> vlan2time;
    vlan2time.clear();
    vector<int> x0;// randomly allocate at the beginning
    x0.clear();
    for(int i=0;i<numhash;i++)
        x0.push_back(i);
//    for(int i=0;i<x0.size();i++)
//        cout<<x0[i]<<" ";
    for(int i=0;i<numhash;i++){
        int t = x0[i];
        x0[i] = t%numvlan ;
    }
    // now x0 = 0,1,2,3,0,1,2,3
    int assign_time=100;
    int t=0;
    //每个vlan assignment的停留时间
    //test
    double cost_old = inf;
    while(true){
        string assign_now = vector2string(x0);
//        cout<<"assign_now = "<<assign_now<<endl;
//        cout<<"x0 = "<<endl;
//        for(int i=0;i<numhash;i++)
//            cout<<x0[i]<<" ";
//        cout<<endl;
        double stay_time_now = 0.0;//time staying in this state(vlan assignment)
        // compute network cost
        MatrixXd R,rate;
        R = getRoute(problem,x0);
        //cout<<"R = "<<endl<<R<<endl;
        MatrixXd active = tm;
        for(int i=0;i<active.rows();i++){
            if(active(i,0)>0)
                active(i,0) = 1;
        }
        rate = maxminRate(problem,active,R);

        MatrixXd links_utilization;
        links_utilization = R*rate;

        double cost_now = 0;
        for(int i=0;i<problem.num_link;i++){
            links_utilization(i,0) = links_utilization(i,0)/problem.capacity(i,0);
//            if(links_utilization(i,0)>cost_now)
//                cost_now = links_utilization(i,0); //find the most congested link

            cost_now += links_utilization(i,0)*links_utilization(i,0);

        }
        cout<<"cost_old = "<<cost_old<<endl;
        cout<<"cost_now = "<<cost_now<<endl;
        if(fabs(cost_now-cost_old)<1e-5)
            break;
        cost_old = cost_now;
        // generate numhash timers(random numbers)
        //vector<double> timers;
        double lambda = alpha*(numvlan-1)*exp(beta*cost_now);
        //cout<<"lambda = "<<lambda<<endl;

        std::exponential_distribution<double> distribution(lambda);

        int index_of_min_timer;//index of the hash class
        double min_timer = 1000000;//stay time
        for(int i=0;i<numhash;i++){
            double r;
            r = distribution(generator);
            //cout<<"r = "<<r<<endl;
            if(r < min_timer){
                min_timer = r;
                index_of_min_timer = i;
               //
            }

        }
        //keep a record of the min_timer
        stay_time_now = min_timer;
        //test
        //cout<<"assign = "<<assign_now<<"  min_timer = "<<min_timer<<endl;
        //cout<<"assign = "<<assign_now<<"  index_of_min_timer = "<<index_of_min_timer<<endl;
        map<string,double>::iterator it = vlan2time.find(assign_now);
        if(it != vlan2time.end()){
            it->second += stay_time_now;
        }else{
            vlan2time[assign_now] = stay_time_now;
        }

        // the hash class : index_of_min_timer needs to be reassigned
        int origin_vlan = x0[index_of_min_timer];
        // randomly pick from the left vlans
        vector<int> left_vlans;
        left_vlans.clear();
        for(int i=0;i<numvlan;i++){
            if(i == origin_vlan)
                continue;
            left_vlans.push_back(i);
        }
        std::uniform_int_distribution<int> unifor_distribution(0,numvlan-2);
        int next_vlan = left_vlans[unifor_distribution(generator)];
        x0[index_of_min_timer] = next_vlan;
        t++;
    }

    // find the vlan assignment with the longest time
    double max_time = -1;
    string vlan="";
    for(map<string,double>::iterator it = vlan2time.begin();it!=vlan2time.end();it++){
        if(it->second > max_time){
            max_time = it->second;
            vlan = it->first;
        }
    }

    vector<int> vlan_assignment_final;// the vlan assignment returned
    vlan_assignment_final.clear();
    //计算vlan是几位数
    for(int i=0;i<vlan.size();i++){
        int now = vlan[i] - 48;
        vlan_assignment_final.push_back(now);
    }
    cout<<"vlan end"<<endl;
    return vlan_assignment_final;
}
// x: vlan_assignment
MatrixXd Run::getRoute(Prob problem,vector<int> x){

    int numlink = problem.num_link;
    int numhash = problem.num_routing_class,numhost = problem.num_edge_switch;
    MatrixXd R(numlink,len);

    for(int h=0;h<numhash;h++){
        int v = x[h];
        for(int i=0;i<numlink;i++){
            for(int j=h*numhost*numhost;j<(h+1)*numhost*numhost;j++){
                R(i,j) = problem.link_flow[v](i,j-h*numhost*numhost);
            }
        }
    }
    return R;
}
MatrixXd Run::maxminRate(Prob problem,MatrixXd x,MatrixXd R){
    MatrixXd C = problem.capacity;
    int L = C.rows();
    //cout<<"L = "<<L<<endl;
    int F = x.rows();
    MatrixXd C_res = C;
    MatrixXd f = x;
    MatrixXd rate;
    rate = MatrixXd::Zero(F,1);
    while(f.sum()>=1){
        MatrixXd fcount;
        fcount = R*f;
       // cout<<"fcount = "<<endl<<fcount<<endl;
        MatrixXd fshare;
        fshare = MatrixXd::Ones(L,1);
        fshare = inf*fshare;
        for(int l=0;l<L;l++){
            if(fcount(l,0)>0){
                fshare(l,0) = C_res(l,0)*1.0/fcount(l,0);
            }
        }
        int minlink,c;
        double val = fshare.minCoeff(&minlink,&c);
        for(int i=0;i<F;i++){
            if(R(minlink,i) == 1 && f(i,0)>0){
                rate(i,0) = val;
                f(i,0) = 0;
            }
        }
        C_res = C - R*rate;
    }
    return rate;
}
string Run::vector2string(vector<int> x){
    string s = "";
    for(int i=0;i<(x.size());i++){
        char t = x[i]+48;
        s+=t;
    }
    return s;
}
double Run::get_total_time(){
    return total_time;
}
vector<double> Run::get_total_thr(){
    return total_thr;
}
vector<double> Run::get_opt_time(){
    return opt_time;
}

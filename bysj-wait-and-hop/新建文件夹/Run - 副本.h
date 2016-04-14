#include<iostream>
#include<Eigen/Dense>
#include<vector>
#include<string>
#include<map>
using namespace std;
using namespace Eigen;
class Run{
private:
    Prob problem;
    double total_time;//总完成时间
    double opt_len;
    double beta;
    double alpha;
    MatrixXd tm;
    vector<int> vlan_assignment;
public:
    Run(Prob p,double opt_length);
    double get_total_time();
    MatrixXd get_routings(Prob p,vector<int> vlan_assignment);//若vlan assign = 0 1 2 3，则返回一个矩阵R=[link_flow[0] link_flow[1] ... link_flow[3]]
    MatrixXd maxminRate(Prob problem,MatrixXd x,MatrixXd R);
    vector<int> markov_vlan_assign(Prob p);
    string vector2string(vector<int> x);
};

Run::Run(Prob p,double opt_length){
    beta = 5.0;
    alpha = 1.0;
    opt_len = opt_length;
    problem = p;
    tm = problem.TM;//得到流量矩阵(每列代表一种routing class)
    int len = problem.num_edge_switch*problem.num_edge_switch*problem.num_routing_class;
    tm.resize(len,1);
    vlan_assignment.clear();
    total_time = 0;
    double last_time = 0;
    vlan_assignment = markov_vlan_assign(problem);
    for(int i=0;i<vlan_assignment.size();i++){
        cout<<vlan_assignment[i]<<" ";
    }
    cout<<endl;
    while(1){
        if(opt_len<0){
            vlan_assignment = markov_vlan_assign(problem);

        }else if(total_time>opt_len + last_time){
            vlan_assignment = markov_vlan_assign(problem);
            last_time = total_time;
        }
        MatrixXd active = tm;
        for(int i=0;i<active.rows();i++){
            if(active(i,0)>0)
                active(i,0) = 1;
        }

        MatrixXd R;
        R = get_routings(problem,vlan_assignment);

        MatrixXd rate;
        rate = maxminRate(problem,active,R);

        cout<<active.sum()<<" flows left\n";
        double duration = inf*1.0;
        bool changed = false;
        for(int i=0;i<len;i++){
            if(tm(i,0) > 0){
                double t = tm(i,0)/rate(i,0);
                if(t < duration){
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

    }

}
double Run::get_total_time(){
    return total_time;
}
MatrixXd Run::get_routings(Prob p,vector<int> vlan_assignment){
    MatrixXd r;
    r = MatrixXd::Zero(p.num_link,p.num_flow*p.num_routing_class);
    for(int i=0;i<vlan_assignment.size();i++){
        int vlan = vlan_assignment[i];
        for(int row=0;row<p.num_link;row++){
            for(int col = i*p.num_flow;col<(i+1)*p.num_flow;col++){
                r(row,col) = p.link_flow[vlan](row,col-i*p.num_flow);
            }
        }
    }
    return r;
}

vector<int> Run::markov_vlan_assign(Prob p){
    int num_routing_class = p.num_routing_class;
    int num_vlan = p.num_vlan;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    map<string,double> vlan2time;
    vlan2time.clear();
    vector<int> x0;// randomly allocate at the beginning
    x0.clear();
    for(int i=0;i<num_routing_class;i++)
        x0.push_back(i);
//    for(int i=0;i<x0.size();i++)
//        cout<<x0[i]<<" ";
    for(int i=0;i<num_routing_class;i++){
        int t = x0[i];
        x0[i] = t%num_vlan ;
    }
    // now x0 = 0,1,2,3,0,1,2,3
    int t=0;
    //每个vlan assignment的停留时间
    //test
    double cost_old = 1000000;
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
        R = get_routings(p,x0);
        MatrixXd active = tm;
        for(int i=0;i<active.rows();i++){
            if(active(i,0)>0)
                active(i,0) = 1;
            cout<<active(i,0)<<endl;
        }
        rate = maxminRate(p,active,R);

        MatrixXd links_utilization;
        links_utilization = R*rate;

        double cost_now=0;
        for(int i=0;i<links_utilization.rows();i++){
            links_utilization(i,0) = links_utilization(i,0)/problem.capacity(i,0);
            cost_now += links_utilization(i,0)*links_utilization(i,0);

        }
        if(fabs(cost_now - cost_old)<1e-5)
            break;
        // generate numhash timers(random numbers)
        //vector<double> timers;
        double lambda = alpha*(num_vlan-1)*exp(beta*cost_now);
        cout<<"lambda = "<<lambda<<endl;
        std::exponential_distribution<double> distribution(lambda);

        int index_of_min_timer;//index of the hash class
        double min_timer = 1000000;//stay time
        for(int i=0;i<num_routing_class;i++){
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
        for(int i=0;i<num_vlan;i++){
            if(i == origin_vlan)
                continue;
            left_vlans.push_back(i);
        }
        std::uniform_int_distribution<int> unifor_distribution(0,num_vlan-2);
        int next_vlan = left_vlans[unifor_distribution(generator)];
        x0[index_of_min_timer] = next_vlan;
        t++;
        cout<<"vlan assignment : t = "<<t<<endl;
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
//    cout<<"string vlan = "<<vlan<<endl;
//    cout<<"vlan_assignment = ";
//    for(int i=0;i<vlan_assignment_final.size();i++)
//        cout<<vlan_assignment_final[i]<<" ";
//    cout<<endl;
    //cout<<"hi "<<vlan_assignment_final[0]<<endl;
    return vlan_assignment_final;
}

MatrixXd Run::maxminRate(Prob problem,MatrixXd x,MatrixXd R){
    MatrixXd C = problem.capacity;
    int L = C.rows();
    int F = x.rows();
    MatrixXd C_res = C;
    MatrixXd f = x;
    MatrixXd rate;
    rate = MatrixXd::Zero(F,1);
    while(f.sum()>=1){
        MatrixXd fcount;
        fcount = R*f;
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

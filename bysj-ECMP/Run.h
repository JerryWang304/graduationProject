/*
*  µœ÷ECMPÀ„∑®
*/
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
    int len;
    double total_time;
    double last_time;
    Prob problem;
    vector<int> vlan_assignment;
    MatrixXd tm;
public:
    Run(Prob prob);
    MatrixXd getRoute(Prob problem);
    MatrixXd maxminRate(Prob problem,MatrixXd x,MatrixXd R);
    double get_total_time();

};
Run::Run(Prob prob){
    problem = prob;
    len = problem.num_routing_class*problem.num_edge_switch*problem.num_edge_switch;
    total_time = 0;
    last_time = 0;
    tm = problem.TM;
    tm.resize(len,1);
//    for(int i=0;i<vlan_assignment.size();i++)
//        cout<<vlan_assignment[i]<<" ";
//    cout<<endl;
    while(1){

        MatrixXd active = tm;
        for(int i=0;i<active.rows();i++){
            if(active(i,0)>0)
                active(i,0) = 1;
        }

        MatrixXd R;
        R = getRoute(problem);

        MatrixXd rate;
        rate = maxminRate(problem,active,R);

        cout<<active.sum()<<" flows left\n";
        double duration = inf*1.0;
        bool changed = false;
        for(int i=0;i<len;i++){
            if(tm(i,0) > 0){
                double t = tm(i,0)/rate(i,0);

                if(0<t && t < duration){
                    //cout<<tm(i,0)<<" "<<rate(i,0)<<endl;
                    duration = t;
                    changed = true;
                }
            }
        }
        //cout<<"duration = "<<duration<<endl;
        if(changed == false)
            break;
        total_time += duration;
        for(int i=0;i<len;i++){
            if(tm(i,0) > 0)
                tm(i,0) = tm(i,0) - rate(i,0)*duration;
        }
    }
}

MatrixXd Run::getRoute(Prob problem){

    return problem.link_flow;
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
double Run::get_total_time(){
    return total_time;
}


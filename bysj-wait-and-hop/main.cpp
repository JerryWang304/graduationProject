/**
* 功能：在fat-tree中实现"wait and hop"算法
* 作者：王丰
* 日期：2016-4-1
*/
#include <iostream>
#include"Prob.h"
#include"Run.h"
#include<fstream>
using namespace std;
int main()
{


    ofstream out;
//
//    for(int i=0;i<100;i++){
//        out.open("F:\\wait_and_hop_k_8.txt",ios::app);
//        Prob p(8);
//        Run r(p,1);
//        if(i<99)
//            out<<r.get_total_time()<<" ";
//        else
//            out<<r.get_total_time();
//        cout<<"第"<<i<<"次： "<<r.get_total_time()<<endl;
//        out.close();
//
//    }
    Prob p(4);
    Run r(p,1);
    cout<<r.get_total_time();
    return 0;
}

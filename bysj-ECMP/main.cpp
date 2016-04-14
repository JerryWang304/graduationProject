/**
* 功能：在fat-tree中实现ECMP算法
* 作者：王丰
* 日期：2016-4-5
*/
#include <iostream>
#include"Prob.h"
#include"Run.h"
#include<fstream>
using namespace std;
int main()
{
    ofstream out;
    for(int i=0; i<100; i++)
    {
        out.open("F:\\ECMP-8.txt",ios::app);
        Prob p(8);
        Run r(p);
        if(i<99)
            out<<r.get_total_time()<<" ";
        else
            out<<r.get_total_time();
        out.close();
        cout<<"i = "<<i<<endl;
    }
    return 0;
}

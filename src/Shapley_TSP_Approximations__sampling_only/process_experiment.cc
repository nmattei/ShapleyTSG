#include<cstdlib>
#include<cassert>
#include<sys/types.h>
#include<time.h>

#include<fstream>
#include<map>
#include<vector>
#include<string>
#include<iostream>
#include<sstream>


using namespace std;

vector<double> tmst_small;
vector<double> mst_small;
vector<double> tsp_small;

vector<double> tmst_large;
vector<double> mst_large;
vector<double> tsp_large;


vector<double> mst_average;
vector<double> tmst_average;
vector<double> tsp_average;
vector<double> tsp_mst_average_difference;

double calculate_average_difference(vector<double>& mst_large, vector<double>& mst_small, vector<double>& tsp_large, vector<double>& tsp_small){
    assert(mst_large.size() == mst_small.size());
    assert(tsp_large.size() == tsp_small.size());
    assert(tsp_large.size() == mst_small.size());
    
    double count=0;
    double sum;
    for (size_t i =0 ; i < mst_large.size(); i++){
        double token_mst = mst_large[i] - mst_small[i];
        double token_tsp = tsp_large[i] - tsp_small[i];
        
        // tsp * x = mst :: x = mst / tsp

        //if(token_mst < 0){cerr<<"Warning margin is :: "<<token_mst<<endl;continue;}
        //if(token_tsp < 0){cerr<<"Warning margin is :: "<<token_tsp<<endl;continue;}

        if(token_tsp == 0){cerr<<"."; char ch;cin>>ch;continue;}
        assert(token_tsp != 0);

        //assert(token >= 0);
        //assert(token_mst > token_tsp);
        count+=1;
        
        sum += token_mst / token_tsp;
        //sum += token_mst - token_tsp;//token;
    }
    cerr<<sum<<" "<<count<<endl;

    sum = sum / count;//static_cast<double>(large.size());
    
    return sum;
}

double calculate_average(vector<double>& large, vector<double>& small, bool with_warning = false){
    assert(large.size() == small.size());
    
    double count=0;
    double sum;
    for (size_t i =0 ; i < large.size(); i++){
        double token = large[i] - small[i];
        
        if(token < 0 && with_warning){cerr<<"Warning margin is :: "<<token<<endl;}//continue;}

        //assert(token >= 0);
        count+=1;
        sum += token;
    }
    sum = sum / count;//static_cast<double>(large.size());
    
    return sum;
}

void parse_small_points(ifstream& in){

    if(in.eof())return;
    string tmp;
    do{
        getline(in, tmp);
        if(!tmp.size())continue;
        istringstream iss(tmp);
        string str;
        iss>>str;
        double value;
        if(str == "MST"){
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>value;
            mst_small.push_back(value);
        }
        
        if (str == "TSP"){
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>value;
            tsp_small.push_back(value);
        }
        
        if (str == "TMST"){
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>value;
            tmst_small.push_back(value);
        }
        
    }while(!in.eof());
}

void parse_large_points(ifstream& in){

    if(in.eof())return;
    string tmp;
    do{
        getline(in, tmp);
        if(!tmp.size())continue;
        istringstream iss(tmp);
        string str;
        iss>>str;
        double value;
        if(str == "MST"){
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>value;
            mst_large.push_back(value);
        }
        
        if (str == "TSP"){
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>value;
            tsp_large.push_back(value);
        }

        if (str == "TMST"){
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>str;
            iss>>value;
            tmst_large.push_back(value);
        }
        
    }while(!in.eof());
}

void parse(size_t i){
    {
        ostringstream oss;
        oss<<"output."<<i;
        string filename = oss.str();
        ifstream file;
        file.open (filename.c_str());
        parse_small_points(file);
        file.close();
    }
    {
        ostringstream oss;
        oss<<"output."<<i<<".1";
        string filename = oss.str();
        ifstream file;
        file.open (filename.c_str());
        parse_large_points(file);
        file.close();
    }

    cerr<<"MST margin first. \n";
    double average_mst = calculate_average(mst_large, mst_small);
    
    cerr<<"TMST margin first. \n";
    double average_tmst = calculate_average(tmst_large, tmst_small);
    
    cerr<<"TSP margin second. \n";
    double average_tsp = calculate_average(tsp_large, tsp_small);

    double average_difference = calculate_average_difference(tmst_large, tmst_small, tsp_large, tsp_small);

    mst_average.push_back(average_mst);
    tmst_average.push_back(average_tmst);
    tsp_average.push_back(average_tsp);
    tsp_mst_average_difference.push_back(average_difference);
}

int main(int argc, char** argv){
    
    size_t lower = atoi(argv[1]); 
    size_t upper = atoi(argv[2]); 
    
    for(int i = lower; i <= upper; i++){
        parse(i);
    }
    
    // For octave
    cout<<"x = [ ";
    for(int i = lower; i <= upper; i++){
        cout<<i<<" ";
    }
    cout<<" ] ; "<<endl;
    
    cout<<"mst = [ ";
    for ( size_t i = 0 ; i < mst_average.size(); i++){
        cout<<mst_average[i]<<" ";
    }
    cout<<" ] ; "<<endl;
    
    cout<<"tmst = [ ";
    for ( size_t i = 0 ; i < tmst_average.size(); i++){
        cout<<tmst_average[i]<<" ";
    }
    cout<<" ] ; "<<endl;
    
    cout<<"tsp = [ ";
    for ( size_t i = 0 ; i < tsp_average.size(); i++){
        cout<<tsp_average[i]<<" ";
    }
    cout<<" ] ; "<<endl;
    
    cout<<"tsp_mst = [ ";
    for ( size_t i = 0 ; i < tsp_mst_average_difference.size(); i++){
        cout<<tsp_mst_average_difference[i]<<" ";
    }
    cout<<" ] ; "<<endl;
    
    
    return 0;
}

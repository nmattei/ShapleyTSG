

/*Map from concorde points back to experiment points.*/
map<int, int> map_from_concorde_points_to_local_points;


void write__for_concorde__matrix(ofstream& file, const set<int>& local_points){
    
    map_from_concorde_points_to_local_points = map<int, int>();
    assert(local_points.size() > 3);
    

    file<<"NAME:  br17"<<endl;
    file<<"TYPE: TSP"<<endl;
    file<<"COMMENT: 17 city problem (Repetto)"<<endl;
    file<<"DIMENSION:  "<<local_points.size()<<endl;
    file<<"EDGE_WEIGHT_TYPE: EXPLICIT"<<endl;
    file<<"EDGE_WEIGHT_FORMAT: FULL_MATRIX"<<endl;
    file<<"EDGE_WEIGHT_SECTION"<<endl;
    
    int i = 1;
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        
        //cerr<<i-1<<" "<<*p<<endl;
        map_from_concorde_points_to_local_points[i-1] = *p;
        i++;
        
        for(set<int>::const_iterator q = local_points.begin()
            ; q != local_points.end()
            ; q++){
            //cerr<<COST_OF_EDGE__FUNCTION(*p, *q)<<" ";//<<endl;
            //cout<<COST_OF_EDGE__FUNCTION(*p, *q)<<endl;
            // 11 
            file<<(static_cast<int>(COST_OF_EDGE__FUNCTION(*p, *q)))<<" ";//matrix[*p][*q]<<" ";
        }
        file<<endl;

        //cerr<<endl;
    }
    file<<"EOF"<<endl;
    //exit(0);
}


void write__for_concorde__matrix(const set<int>& local_points){
    ofstream file;
    file.open (concorde_file_name.c_str());
    write__for_concorde__matrix(file, local_points);
    file.close();
}

void write__for_concorde(ofstream& file, const set<int>& local_points){
    map_from_concorde_points_to_local_points = map<int, int>();
    assert(local_points.size() > 3);

    file<<"NAME: shapley.tsp"<<endl;
    file<<"TYPE: TSP"<<endl;
    file<<"COMMENT: Shapley Experiment"<<endl;
    file<<"DIMENSION: "<<local_points.size()<<endl;
    file<<"EDGE_WEIGHT_TYPE: EUC_2D"<<endl;
    file<<"DISPLAY_DATA_TYPE: COORD_DISPLAY"<<endl;
    file<<"NODE_COORD_SECTION"<<endl;
    int i = 1;
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        //cerr<<"Concorde mapping :: "<<i-1<<" -to-> "<<*p<<endl;
        map_from_concorde_points_to_local_points[i-1] = *p;
        file<<i++<<" "<<points[*p].first<<" "<<points[*p].second<<endl;
    }
    file<<"EOF"<<endl;
}

void write__for_concorde(const set<int>& local_points){

    ofstream file;
    file.open (concorde_file_name.c_str());
    write__for_concorde(file, local_points);
    file.close();
}




void write__for_indigo__matrix(ofstream& file, const set<int>& local_points){
    cerr<<"Writing for indigo. \n";
    
    map_from_concorde_points_to_local_points = map<int, int>();
    assert(local_points.size() > 3);
    
    file<<"VRX"<<endl;
    file<<"NAME Real-world scenario from TipTop2013."<<endl;
    file<<"COMMENT A route of length 20 or 10."<<endl;
    
    
    file<<"COMMODITIES Load"<<endl;
    file<<"METRICS"<<endl;
    file<<"  Time MAT MULT 100"<<endl;
    file<<"*END*"<<endl;
    file<<"TIME_METRIC Time"<<endl;
    
    
    file<<"LOCATIONS"<<endl;
    
    int i = 0;
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        
        //cerr<<i-1<<" "<<*p<<endl;
        map_from_concorde_points_to_local_points[i] = *p;
        
        file<<i<<endl;
        i++;
    }
    file<<"*END*"<<endl;
    
    
    
    file<<"VEHICLES"<<endl;
    file<<"  V1 1000"<<endl;
    file<<"*END*"<<endl;
    
    file<<"VEHICLE_AVAIL"<<endl;
    file<<"  V1 0 0 0 10000"<<endl;
    file<<"*END*"<<endl;
    
    
    file<<"REQUESTS"<<endl;
    i = 0;
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        
        file<<i<<" "<<i<<" 1000 1"<<endl;
        i++;
    }
    file<<"*END*"<<endl;
    
    file<<"EDGE_WEIGHT Time"<<endl;
    i = 0;
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        
        
        int j = 0;
        for(set<int>::const_iterator q = local_points.begin()
            ; q != local_points.end()
            ; q++){

            if(i != j) assert((static_cast<int>(COST_OF_EDGE__FUNCTION(*p, *q)) ) != 0);

            file<<i<<" "<<j<<" "<<(static_cast<int>(COST_OF_EDGE__FUNCTION(*p, *q)) )<<endl;//<<" ";//matrix[*p][*q]<<" ";

            j++;
        }
        
        i++;
        //cerr<<endl;
    }
    file<<"*END*"<<endl;
    
}
    
void write__for_indigo__matrix(const set<int>& local_points){

    ofstream file;
    file.open (concorde_file_name.c_str());
    write__for_indigo__matrix(file, local_points);
    file.close();
}



/* Return value :: First element is total distance, and second element is list of
 * marginal distances for each possible single-customer removal. */
pair<double, vector<double> >  parse_tour(ifstream& in, double(*distance_calc)(int i, int j))
{
    //cerr<<"Parsing tour"<<endl;
    
    if(in.eof()){assert(0);exit(-1);return pair<double, vector<double> >();}
    string tmp;
    getline(in, tmp);

    if(in.eof()){assert(0);exit(-1);return pair<double, vector<double> >();}
    
    double total_distance=0;

    
    size_t first_vertex_id = 1000000;
    size_t old_vertex_id = 1000000;
    size_t _vertex_id;
    size_t vertex_id;
    
    vector<int> route;
    do{
        getline(in, tmp);
        if(in.eof())break;
        if(!tmp.size())continue;
        //cerr<<tmp<<endl;

        istringstream iss(tmp);
        
        do {
            iss>>_vertex_id;
            
            if(map_from_concorde_points_to_local_points.find(_vertex_id) 
               == map_from_concorde_points_to_local_points.end()){
                cerr<<"Error cannot find Concorde vertex :: "<<_vertex_id<<endl;

                //cerr<<
            }

            assert(map_from_concorde_points_to_local_points.find(_vertex_id) 
                   != map_from_concorde_points_to_local_points.end());

            vertex_id =map_from_concorde_points_to_local_points[_vertex_id];


            if(iss.eof())break;
            route.push_back(vertex_id);
            
            if(old_vertex_id == 1000000){
                old_vertex_id = vertex_id;

                first_vertex_id = vertex_id; 
                continue;
            } else {
                //cerr<<old_vertex_id<<" to "<<vertex_id<<endl;

                double distance = distance_calc(old_vertex_id, vertex_id);
                
                total_distance+=distance;
            }

            old_vertex_id = vertex_id;
        }while(!iss.eof());
        
        
    }while(!in.eof());

    double distance = distance_calc(vertex_id, first_vertex_id);
    total_distance+=distance;
    
    double tsp_tour_length = total_distance; 
    
    vector<double> cost_allocations(route.size());
    for(size_t i =0 ; i < route.size(); i++){
        int point_index = route[i];
        
        if(are_we_solving_a_problem_with_a_single_depot){
            if(point_index == 0){
                cost_allocations[point_index] = 0.0;
                continue;
            }
        }
        
        double allocation = 0;
        int previous_point = (i>0)?route[i-1]:route.back();
        int next_point = (i != route.size() -1)?route[i+1]:*route.begin();
        // cerr<<"Distance from : "<<previous_point<<" to "<<point_index
        //     <<" and then from "<<point_index<<" to "<<next_point<<endl;
        double old_distance = distance_calc(previous_point, point_index) + distance_calc(point_index, next_point); 
        double new_distance = distance_calc(previous_point, next_point); 


        //cerr<<old_distance<<" "<<new_distance<<" "<<old_distance - new_distance<<endl;
        assert(i < cost_allocations.size());

        cost_allocations[i] = old_distance - new_distance;

        

        //assert(cost_allocations[i] >= 0);
    }

    //cout<<"TSP Tour length is = "<<total_distance<<endl;
    pair<double, vector<double> > answer(tsp_tour_length, cost_allocations);
    return answer;
}


void parse_Euclidean_tsplib_model(ifstream& in){
    points = vector<Point>();
    
    if(in.eof()){cerr<<"Unrecoverable Error : Empty input file.\n";assert(0);exit(-1);}
    string tmp;
    bool parsing_point_information = false;
    size_t point_counter = 0;
    while(!in.eof()){
        getline(in, tmp);
        istringstream iss(tmp);
        
        string first = "";
        if(!parsing_point_information){
            iss>>first;
        }

        if(!parsing_point_information && "DIMENSION:" == first){
            size_t size;
            iss>>size;
            lower_range = upper_range = points_count = size;
        } else if (!parsing_point_information && "NODE_COORD_SECTION" == first){
            parsing_point_information = true;
        } else if (parsing_point_information){
            double i,j;
            iss>>i>>j;
            points.push_back(Point(i, j));
            point_counter++;
            if(point_counter >= points_count * points_count) break;
        }
    }
}

void parse_matrix_tsplib_model(ifstream& in)
{
    cerr<<"Parsing matrix TSP model.\n";
    if(in.eof()){cerr<<"Unrecoverable Error : Empty input file.\n";assert(0);exit(-1);}
    string tmp;
    bool parsing_point_information = false;
    size_t row_counter = 0;
    while(!in.eof()){
        getline(in, tmp);
        cerr<<tmp<<endl;

        istringstream iss(tmp);
        
        string first = "";
        if(!parsing_point_information){
            iss>>first;
        }
        
        if(!parsing_point_information && "DIMENSION:" == first){
            size_t size;
            iss>>size;
            lower_range = upper_range = points_count = size;
            matrix = vector<vector<double> >(size);

            cerr<<"Problem is of size: "<<points_count<<endl;
        } else if (!parsing_point_information && "EDGE_WEIGHT_SECTION" == first){
            parsing_point_information = true;
        } else if (parsing_point_information){
            assert(points_count>0);
            matrix[row_counter] = vector<double>(points_count);
            size_t entry_counter = 0;
            while(!iss.eof()){
                double item;
                iss>>item;
                //cerr<<item<<" "<<entry_counter<<endl;
                matrix[row_counter][entry_counter++] = item;
                if(entry_counter == points_count) break;
            }

            points.push_back(Point(0,0));

            assert(entry_counter == points_count);
            row_counter++;

            if(row_counter >= points_count) break;
        }
    }
    assert(points.size());
}

void parse_Euclidean_tsplib_model(const string& in){
    ifstream file;
    file.open (in.c_str());
    
    if(file.fail()){
        cerr<<"Unrecoverable Error : Cannot open Euclidean TSP model named "<<in<<endl;
    }

    parse_Euclidean_tsplib_model(file);
    file.close();
}

void parse_matrix_tsplib_model(const string& in){
    ifstream file;
    file.open (in.c_str());
    
    if(file.fail()){
        cerr<<"Unrecoverable Error : Cannot open matrix TSP model named "<<in<<endl;
    }

    parse_matrix_tsplib_model(file);
    file.close();
}

pair<double, vector<double> >  parse_tour__indigo(ifstream& in, double(*distance_calc)(int i, int j))
{
    //cerr<<"Parsing tour"<<endl;
    
    if(in.eof()){assert(0);exit(-1);return pair<double, vector<double> >();}
    string tmp;
    do {
        getline(in, tmp);
    }while(!in.eof() && tmp[0] == '#');

    //We are up-to the "COST" line, now we can parse the route.
    //So the next line is the route.
    getline(in, tmp);

    istringstream iss(tmp);
    string scrap;
    iss>>scrap>>scrap>>scrap;

    double total_distance=0;

    
    size_t first_vertex_id = 1000000;
    size_t old_vertex_id = 1000000;
    size_t _vertex_id;
    size_t vertex_id;
    
    vector<int> route;
    do {
        iss>>_vertex_id;
        
        if(map_from_concorde_points_to_local_points.find(_vertex_id) 
           == map_from_concorde_points_to_local_points.end()){
            cerr<<"Error cannot find Concorde vertex :: "<<_vertex_id<<endl;
            
            //cerr<<
        }
        
        assert(map_from_concorde_points_to_local_points.find(_vertex_id) 
               != map_from_concorde_points_to_local_points.end());
        
        vertex_id =map_from_concorde_points_to_local_points[_vertex_id];
        
        
        if(iss.eof())break;
        route.push_back(vertex_id);
        
        if(old_vertex_id == 1000000){
            old_vertex_id = vertex_id;
            
            first_vertex_id = vertex_id; 
            continue;
        } else {
            //cerr<<old_vertex_id<<" to "<<vertex_id<<endl;
            
            double distance = distance_calc(old_vertex_id, vertex_id);
            
            total_distance+=distance;
        }
        
        old_vertex_id = vertex_id;
    }while(!iss.eof());

    double distance = distance_calc(vertex_id, first_vertex_id);
    total_distance+=distance;
    
    double tsp_tour_length = total_distance; 
    
    vector<double> cost_allocations(route.size());
    for(size_t i =0 ; i < route.size(); i++){
        int point_index = route[i];
        
        if(are_we_solving_a_problem_with_a_single_depot){
            if(point_index == 0){
                cost_allocations[point_index] = 0.0;
                continue;
            }
        }
        
        double allocation = 0;
        int previous_point = (i>0)?route[i-1]:route.back();
        int next_point = (i != route.size() -1)?route[i+1]:*route.begin();
        // cerr<<"Distance from : "<<previous_point<<" to "<<point_index
        //     <<" and then from "<<point_index<<" to "<<next_point<<endl;
        double old_distance = distance_calc(previous_point, point_index) + distance_calc(point_index, next_point); 
        double new_distance = distance_calc(previous_point, next_point); 


        //cerr<<old_distance<<" "<<new_distance<<" "<<old_distance - new_distance<<endl;
        assert(i < cost_allocations.size());

        cost_allocations[i] = old_distance - new_distance;

        

        //assert(cost_allocations[i] >= 0);
    }

    //cout<<"TSP Tour length is = "<<total_distance<<endl;
    pair<double, vector<double> > answer(tsp_tour_length, cost_allocations);
    return answer;
}





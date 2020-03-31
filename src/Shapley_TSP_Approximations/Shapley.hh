
void updateShapley(int elements, 
                   vector<int>& permutation, 
                   vector<mpf_class>& answers, 
                   double(*calculate_value_of)(const set<int> &))
{
    uint j = 0;
    if(are_we_solving_a_problem_with_a_single_depot){
        assert(permutation[0] == 0);
        mpf_class& answer = answers[permutation[0]];
        answer = 0.0;
        j = 1;
    }

    for(; j < elements; j++){
        assert(permutation[j] < answers.size());
        int element_index = permutation[j]; 
        mpf_class& answer = answers[element_index];
            
        set<int> prefix;
            
        int index =0; 
        if(are_we_solving_a_problem_with_a_single_depot){
            prefix.insert(permutation[0]);
            assert(permutation[0] == 0);
        }
        while(permutation[index] != element_index){
            prefix.insert(permutation[index]);
            index++;
        }

        double prefix_cost = 0.0;
        if(prefix.size() > 1){
            //cerr<<"Calculating for prefix : \n";
            prefix_cost = calculate_value_of(prefix);
        }
        
        set<int>& suffix = prefix; // Last customer is added to the prefix to yield a suffix.
        suffix.insert(permutation[index]);
        mpf_class suffix_cost = 0.0;
        
        if(suffix.size() > 1){
            //cerr<<"Calculating for suffix : \n";
            suffix_cost = calculate_value_of(suffix);
        }

        mpf_class delta = (suffix_cost - prefix_cost);
        
        if(tspmode && delta < 0 ) {cerr<<"Prefix size : "<<suffix.size()-1
                                       <<" Suffix size : "<<suffix.size()
                                       <<" "<<delta<<endl;}// assert(delta > -5);delta = 0;}
        if(delta < 0 ) delta = 0; // APRIL 1st
        // assert(delta >= 0);
        answer += delta;
    }
}

vector<mpf_class> ApproShapley(
    int elements,
    int samples_count, 
    double(*calculate_value_of)(const set<int> &))
{
    set<int> all_points;
    vector<int> permutation;
    for(uint i = 0 ; i < elements; i++){
        permutation.push_back(i);
        all_points.insert(i);
    }
    
    assert(permutation[0]==0);
    
    vector<mpf_class> answers(elements);
    for(uint sample_number = 0; sample_number < samples_count; sample_number++){
        if(are_we_solving_a_problem_with_a_single_depot){
            vector<int>::iterator begin_at = permutation.begin();
            begin_at++;
            std::random_shuffle(begin_at, permutation.end());
            assert(permutation[0]==0);
        } else {
            std::random_shuffle(permutation.begin(), permutation.end());
        }
        updateShapley(elements, permutation, answers, calculate_value_of);
    }
    mpf_class v_value = calculate_value_of(all_points);//whole_cost;
    mpf_class u_value = 0;
    for(uint j = 0; j < elements; j++){
        u_value +=  answers[j];
    }

    assert(u_value != 0);
    
    for(uint j = 0; j < elements; j++){
        mpf_class& answer = answers[j];
        //cerr<<answer<<" * v/u "<<u_value<<" "<<v_value<<std::endl;
        answer = answer * (v_value / u_value);
        cerr<<"Got a shapley value : "<<answer<<std::endl;
    }

    return answers;
}

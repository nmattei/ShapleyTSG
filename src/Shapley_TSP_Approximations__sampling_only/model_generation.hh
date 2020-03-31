
void generate_model_matrix(){
    matrix = vector<vector<double> >(points_count);
    for (size_t i = 0 ; i < points_count; i++){
        matrix[i] = vector<double>();
        for (size_t j = 0 ; j < points_count; j++){
            if(i == j){
                matrix[i].push_back(ceil(999.0));
            } else {
                if( j < i){
                    if(random() % 2){
                        matrix[i].push_back(ceil( matrix[j][i] * 1.5));
                    } else {
                        matrix[i].push_back(ceil(matrix[j][i] * 0.5));
                    }
                } else {
                    matrix[i].push_back(ceil(drand48() * chart_magnitude));
                }
            }
        }
    }
}

void generate_model(){
    points = vector<Point>();
    double x;
    double y;
    
    for (size_t i =0 ; i < points_count; i++){
        x = drand48() * chart_magnitude;
        y = drand48() * chart_magnitude;
        points.push_back(Point(x, y));
    }
}

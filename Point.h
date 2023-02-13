#ifndef K_MEANS_MIO_CPP_POINT_H
#define K_MEANS_MIO_CPP_POINT_H

class Point {

public:
    Point(double x_coord){
        this->x_coord = x_coord;
        cluster_id = 0;
    }

    Point(){
        x_coord = 0;
        cluster_id = 0;
    }

    double set_x_coord(double x_coord){
        this->x_coord = x_coord;
        cluster_id = 0;
    }

    double get_x_coord(){
        return this->x_coord;
    }


    int get_cluster_id(){
        return cluster_id;
    }

    void set_cluster_id(int cluster_id){
        this->cluster_id = cluster_id;
    }

private:
    double x_coord;
    int cluster_id;
};


#endif //K_MEANS_MIO_CPP_POINT_H

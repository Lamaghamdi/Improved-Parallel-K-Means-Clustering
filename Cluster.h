#ifndef K_MEANS_MIO_CPP_CLUSTER_H
#define K_MEANS_MIO_CPP_CLUSTER_H


#include <queue>
#include "Point.h"
#include <omp.h>

class Cluster
{
public:
    Cluster(double x_coord)
    {
        new_x_coord = 0;
        size = 0;
        this->x_coord = x_coord;
    }

    Cluster()
    {
        new_x_coord = 0;
        size = 0;
        this->x_coord = 0;
    }


    void add_point(Point point)
    {
        #pragma omp atomic
        new_x_coord += point.get_x_coord();
        #pragma omp atomic
        size++;
    }

    void free_point()
    {
        this->size = 0;
        this->new_x_coord = 0;
    }

    double get_x_coord()
    {
        return this->x_coord;
    }

    void set_x_coord(double x_coord)
    {
        this->x_coord = x_coord;
    }

    int get_size()
    {
        return this->size;
    }


    bool update_coords()
    {

        if(this->x_coord == new_x_coord/this->size)
        {
            double x =new_x_coord/this->size;
            return false;
        }

        this->x_coord = new_x_coord/this->size;

        return true;

    }

private:
    double x_coord;
    double new_x_coord;
    int size;
};


#endif //K_MEANS_MIO_CPP_CLUSTER_H

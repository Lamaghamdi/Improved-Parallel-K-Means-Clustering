#include <iostream>
#include "math.h"
#include <fstream>
#include <chrono>
#include "Point.h"
#include "Cluster.h"
#include "stdio.h"
#include <omp.h>
#include<windows.h>

using namespace std;
using namespace std::chrono;

double max_range = 100000;
int num_point = 200000;
int num_cluster =  20;
int max_iterations = 20;


vector<Point> init_point(int num_point);
vector<Cluster> init_cluster(int num_cluster);

void add_point_cluster(int num_cluster, int num_point,vector<Point> &points,vector<Cluster> &clusters);
void set_centroids(int num_cluster, int num_point,vector<Point> &points,vector<Cluster> &clusters);

void compute_distance_euclidean(vector<Point> &points, vector<Cluster> &clusters);
void compute_distance_manhattan(vector<Point> &points, vector<Cluster> &clusters);

void mean_sd(vector<Point> &points, int num_point, double *point_mean,double *point_sd);
double compute_q(vector<Point> &points, int num_point, double *point_mean,double *point_sd);


double euclidean_dist(Point point, Cluster cluster);
double manhattan_dist(Point point, Cluster cluster);

bool update_clusters(vector<Cluster> &clusters);
void draw_chart_gnu(vector<Point> &points);

void print_point_clusters(int num_cluster, int num_point,vector<Point> &points,vector<Cluster> &clusters);

int main()
{

    printf("Number of points %d\n", num_point);
    printf("Number of clusters %d\n", num_cluster);
    printf("Number of processors: %d\n", omp_get_num_procs());

    srand(int(time(NULL)));

    double time_point1 = omp_get_wtime();
    printf("Starting initialization..\n");

    vector<Point> points;
    vector<Cluster> clusters;

    //The Cluster Initialization and the Point initialization are in parallel
    printf("\nCreating points..\n");
    points = init_point(num_point);
    printf("Points initialized \n");

    printf("\nCreating clusters..\n");
    clusters = init_cluster(num_cluster);
    printf("Clusters initialized \n");

    double time_point2 = omp_get_wtime();
    double duration = time_point2 - time_point1;
    printf("\nPoints and clusters generated in: %f seconds\n", duration);

    add_point_cluster(num_cluster,num_point,points,clusters);

    printf("\n***************** Set Initial Centroids *****************\n");
    set_centroids(num_cluster,num_point,points,clusters);

    bool conv = true;
    int iterations = 0;

    double mean, sd;
    mean_sd(points, num_point, &mean, &sd);

    double q = compute_q(points,  num_point, &mean, &sd);

    printf("Starting iterate...\n");

    if (q <= 8.4595)
    {
        printf("\n***************** Compute Distance Euclidean *****************\n");
    }
    else
    {
        printf("\n***************** Compute Distance Manhattan *****************\n");
    }


    //The algorithm stops when iterations > max_iteration or when the clusters didn't move  && iterations < max_iterations
    while(conv && iterations < max_iterations)
    {
        iterations ++;

        if (q <= 8.4595)
        {
            compute_distance_euclidean(points, clusters);
            set_centroids(num_cluster,num_point,points,clusters);
        }
        else
        {
            compute_distance_manhattan(points, clusters);
            set_centroids(num_cluster,num_point,points,clusters);
        }



        conv = update_clusters(clusters);

        printf("Iteration %d done \n", iterations);
    }

    double time_point3 = omp_get_wtime();
    duration = time_point3 - time_point1;

    printf("Number of iterations: %d, total time: %f seconds, time per iteration: %f seconds\n",
           iterations, duration, duration/iterations);

    print_point_clusters(num_cluster,num_point,points,clusters);

    try
    {
        printf("\n\nDrawing the chart...\n\n");
        draw_chart_gnu(points);
    }
    catch(int e)
    {
        printf("Chart not available, gnuplot not found");
    }

    return 0;
}

//Initialize num_point Points
vector<Point> init_point(int num_point)
{
    printf("\n***************** Point initialization *****************\n");

    vector<Point> points(num_point);
    Point *ptr = &points[0];

    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < num_point; i++)
        {
            Point* point = new Point();

            ptr[i] = *point;
        }
    }

    FILE *file = fopen("Employee_Payroll.txt", "r");

    if (file != NULL)
    {
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        int i=0;
        double point_file;


        while ((read = getline(&line, &len, file)) != -1 && num_point>i)
        {
            sscanf(line, "%lf", &point_file);

            points[i].set_x_coord(point_file);
            i++;
        }

        fclose(file);
    }
    else
        printf("Error in opening file\n");

    return points;
}


//Initialize num_cluster Clusters
vector<Cluster> init_cluster(int num_cluster)
{
    printf("\n***************** Cluster initialization *****************\n");

    vector<Cluster> clusters(num_cluster);
    Cluster* ptr = &clusters[0];

    if(num_cluster < omp_get_num_procs())
    {
        #pragma omp parallel num_threads(num_cluster)
        {
            #pragma omp for
            for(int i = 0; i < num_cluster; i++)
            {

                Cluster *cluster = new Cluster();

                ptr[i] = *cluster;
            }
        }
    }
    else
    {
        #pragma omp parallel
        {
            #pragma omp for
            for(int i = 0; i < num_cluster; i++)
            {

                Cluster *cluster = new Cluster();

                ptr[i] = *cluster;

                printf("cluster %d created \n", i+1);

            }
        }
    }


    return clusters;
}


//Initialize num_cluster Clusters
void add_point_cluster(int num_cluster, int num_point, vector<Point> &points, vector<Cluster> &clusters)
{
    printf("\n***************** Add Points *****************\n");

    int initial_clusters_size = floor(num_point/num_cluster);

    int forloop = initial_clusters_size;

    int j = 0;

    for(int i = 0; i < num_cluster; i++)
    {
        Cluster &cluster = clusters[i];

        for (j; j < forloop ; j++)
        {
            Point &point = points[j];
            points[j].set_cluster_id(i);
            clusters[i].add_point(points[j]);
        }

        j=forloop;
        forloop=forloop+initial_clusters_size;
    }
}


//set initial centroids
void set_centroids(int num_cluster, int num_point,vector<Point> &points,vector<Cluster> &clusters)
{
    printf("\n***************** Set Centroids *****************\n");

    for(int i = 0; i < num_cluster; i++)
    {
        Cluster &cluster = clusters[i];
        double point_sum =0;
        int point_counter =0;
        double cluster_mean =0;

        for (int j = 0; j < num_point ; j++)
        {
            Point &point = points[j];

            if (point.get_cluster_id() == i)
            {
                point_sum += point.get_x_coord();
                point_counter += 1;
            }
        }

        cluster_mean = point_sum/point_counter;
        clusters[i].set_x_coord(cluster_mean);
    }
}


//For each Point, compute the distance between each Cluster and assign the Point to the nearest Cluster
//The distance is compute through Euclidean Distance
//The outer for is parallel, with private=min_distance, min_index, points_size, clusters_size and clustes while the
//vector of Points is shared. The amount of computation performed per Point is constant, so static thread scheduling was chosen

void compute_distance_euclidean(vector<Point> &points, vector<Cluster> &clusters)
{
    //printf("\n***************** Compute Distance Euclidean *****************\n");

    unsigned long points_size = points.size();
    unsigned long clusters_size = clusters.size();

    double min_distance;
    int min_index;

    #pragma omp parallel default(shared) private(min_distance, min_index) firstprivate(points_size, clusters_size)
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < points_size; i++)
        {
            Point &point = points[i];

            min_distance = euclidean_dist(point, clusters[0]);
            min_index = 0;

            for (int j = 1; j < clusters_size; j++)
            {

                Cluster &cluster = clusters[j];

                double distance = euclidean_dist(point, cluster);

                if (distance < min_distance)
                {

                    min_distance = distance;
                    min_index = j;
                }

            }
            point.set_cluster_id(min_index);
            clusters[min_index].add_point(point);
        }
    }
}

void compute_distance_manhattan(vector<Point> &points, vector<Cluster> &clusters)
{
    //printf("\n***************** Compute Distance Manhattan *****************\n");

    unsigned long points_size = points.size();
    unsigned long clusters_size = clusters.size();

    double min_distance;
    int min_index;

    #pragma omp parallel default(shared) private(min_distance, min_index) firstprivate(points_size, clusters_size)
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < points_size; i++)
        {
            Point &point = points[i];

            min_distance = manhattan_dist(point, clusters[0]);
            min_index = 0;

            for (int j = 1; j < clusters_size; j++)
            {

                Cluster &cluster = clusters[j];

                double distance = manhattan_dist(point, cluster);

                if (distance < min_distance)
                {

                    min_distance = distance;
                    min_index = j;
                }

            }
            point.set_cluster_id(min_index);
            clusters[min_index].add_point(point);
        }
    }
}

void mean_sd(vector<Point> &points, int num_point, double *point_mean, double *point_sd)
{
    printf("\n***************** Compute Mean and SD *****************\n");

    double point_sum =0;
    int point_counter =0;
    *point_mean =0;
    *point_sd =0;

    for (int i = 0; i < num_point ; i++)
    {
        Point &point = points[i];

        point_sum += point.get_x_coord();
        point_counter += 1;

    }

    *point_mean = point_sum/point_counter;

    for (int i = 0; i < point_counter; ++i)
    {
        Point &point = points[i];
        *point_sd += pow(point.get_x_coord() - *point_mean, 2);
    }

    *point_sd = sqrt(*point_sd / point_counter);
}

double compute_q(vector<Point> &points, int num_point, double *point_mean,double *point_sd)
{

    double q = 0;
    double numerator = 0;
    double denominator = 0;
    int point_counter =0;

    for (int i = 0; i < num_point ; i++)
    {
        Point &point = points[i];

        numerator += pow((point.get_x_coord()- *point_mean), 4);
        point_counter += 1;
    }

    numerator = numerator/point_counter;
    denominator = pow(*point_sd,4);

    q = numerator/denominator;

    return q;
}

double euclidean_dist(Point point, Cluster cluster)
{

    double distance = sqrt(pow(point.get_x_coord() - cluster.get_x_coord(),2));

    return distance;
}

double manhattan_dist(Point point, Cluster cluster)
{

    double distance = abs(point.get_x_coord() - cluster.get_x_coord());

    return distance;
}

//For each cluster, update the coords. If only one cluster moves, conv will be TRUE
//A parallel for was chosen for each cluster with lastprivate=conv
bool update_clusters(vector<Cluster> &clusters)
{

    bool conv = false;

    for(int i = 0; i < clusters.size(); i++)
    {
        conv = clusters[i].update_coords();
        clusters[i].free_point();
    }

    return conv;
}


//Draw point plot with gnuplot
void draw_chart_gnu(vector<Point> &points)
{

    ofstream outfile("data.txt");

    for(int i = 0; i < points.size(); i++)
    {

        Point point = points[i];
        outfile << i << " " << point.get_x_coord() << " " << point.get_cluster_id() << std::endl;

    }

    outfile.close();
    system("gnuplot -p -e \"plot 'data.txt' using 1:2:3 with points palette notitle \"");
    remove("data.txt");

}


void print_point_clusters(int num_cluster, int num_point,vector<Point> &points,vector<Cluster> &clusters)
{

    printf("\n***************** Points In Each Cluster *****************\n");

    vector<int> NoPC(num_cluster);
    double min_max[2];
    unsigned long points_size = points.size();
    unsigned long clusters_size = clusters.size();

    for (int j = 0; j < clusters_size; j++)
    {
        Cluster &cluster = clusters[j];
        min_max[0]= {0.0};
        min_max[1]= {0.0};

        printf("\n************ Cluster %d ************\n", j+1 );

        for (int i = 0; i < points_size ; i++)
        {
            Point &point = points[i];


            if (points[i].get_cluster_id() == j)
            {
                NoPC[j]+=1;

                if (NoPC[j]==1)
                {
                    min_max[0]= points[i].get_x_coord();
                    min_max[1]= points[i].get_x_coord();
                }


                //minimum value in a cluster
                if(points[i].get_x_coord()< min_max[0])
                {
                    min_max[0] = points[i].get_x_coord();
                }


                //maximum value in a cluster
                if(points[i].get_x_coord()> min_max[1])
                {
                    min_max[1] = points[i].get_x_coord();
                }
            }

        }
        printf("\nThere are %d points in cluster %d\n",  NoPC[j], j+1);
        printf("Cluster %d Min is %f Max is %f\n\n", j+1, min_max[0], min_max[1]);

    }
}

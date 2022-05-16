#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 2
#define T0 50000.0  // initial temperature
#define T_end (1e-8)
#define q 0.98      // SA coefficient
#define L 1000      // The number of iterations at each temperature

//#define SELECT_NUM 32       // number of cities is 32
#define SELECT_NUM 64       // number of cities is 64
//#define SELECT_NUM 128      // number of cities is 128
//#define SELECT_NUM 256      // number of cities is 256

#define N SELECT_NUM        // number of city

//global city list consists of 32 different city
//double city_pos[N][2] =
//        {
//                {1304, 2312},
//                {3639, 1315},
//                {4177, 2244},
//                {3712, 1399},
//                {3488, 1535},
//                {3326, 1556},
//                {3238, 1229},
//                {4196, 1004},
//                {4312, 7901},
//                {4386, 5702},
//                {3007, 1970},
//                {2562, 1756},
//                {2788, 1491},
//                {2381, 1676},
//                {1332, 6953},
//                {3715, 1678},
//                {3918, 2179},
//                {4061, 2370},
//                {3780, 2212},
//                {3676, 2578},
//                {4029, 2838},
//                {4263, 2931},
//                {3429, 1908},
//                {3507, 2367},
//                {3394, 2643},
//                {3439, 3201},
//                {2935, 3240},
//                {3140, 3550},
//                {2545, 2357},
//                {2778, 2826},
//                {2370, 2975},
//                {4520, 3412}
//        };

double city_pos[N][2] =
{
		{1304, 2312},
		{3639, 1315},
		{4177, 2244},
		{3712, 1399},
		{3488, 1535},
		{3326, 1556},
		{3238, 1229},
		{4196, 1004},
		{4312, 7901},
		{4386, 5702},
		{3007, 1970},
		{2562, 1756},
		{2788, 1491},
		{2381, 1676},
		{1332, 6953},
		{3715, 1678},
		{3918, 2179},
		{4061, 2370},
		{3780, 2212},
		{3676, 2578},
		{4029, 2838},
		{4263, 2931},
		{3429, 1908},
		{3507, 2367},
		{3394, 2643},
		{3439, 3201},
		{2935, 3240},
		{3140, 3550},
		{2545, 2357},
		{2778, 2826},
		{2370, 2975},
		{4520, 3412},

		{1374, 2212},
		{3679, 1215},
		{4197, 2444},
		{3772, 1299},
		{3478, 1235},
		{3376, 1256},
		{3278, 1429},
		{4176, 1204},
		{4372, 7201},
		{4376, 5202},
		{3077, 1270},
		{2572, 1256},
		{2778, 1291},
		{2371, 1276},
		{1372, 6253},
		{3775, 1278},
		{3978, 2279},
		{4071, 2270},
		{3770, 2412},
		{3696, 2278},
		{4079, 2238},
		{4273, 2231},
		{3449, 1208},
		{3557, 2267},
		{3374, 2243},
		{3419, 3401},
		{2955, 3440},
		{3150, 3250},
		{2575, 2257},
		{2788, 2226},
		{2310, 2275},
		{4530, 3212}
};

//double city_pos[N][2] =
//        {
//                {1304, 2312},
//                {3639, 1315},
//                {4177, 2244},
//                {3712, 1399},
//                {3488, 1535},
//                {3326, 1556},
//                {3238, 1229},
//                {4196, 1004},
//                {4312, 7901},
//                {4386, 5702},
//                {3007, 1970},
//                {2562, 1756},
//                {2788, 1491},
//                {2381, 1676},
//                {1332, 6953},
//                {3715, 1678},
//                {3918, 2179},
//                {4061, 2370},
//                {3780, 2212},
//                {3676, 2578},
//                {4029, 2838},
//                {4263, 2931},
//                {3429, 1908},
//                {3507, 2367},
//                {3394, 2643},
//                {3439, 3201},
//                {2935, 3240},
//                {3140, 3550},
//                {2545, 2357},
//                {2778, 2826},
//                {2370, 2975},
//                {4520, 3412},
//
//                {1374, 2212},
//                {3679, 1215},
//                {4197, 2444},
//                {3772, 1299},
//                {3478, 1235},
//                {3376, 1256},
//                {3278, 1429},
//                {4176, 1204},
//                {4372, 7201},
//                {4376, 5202},
//                {3077, 1270},
//                {2572, 1256},
//                {2778, 1291},
//                {2371, 1276},
//                {1372, 6253},
//                {3775, 1278},
//                {3978, 2279},
//                {4071, 2270},
//                {3770, 2412},
//                {3696, 2278},
//                {4079, 2238},
//                {4273, 2231},
//                {3449, 1208},
//                {3557, 2267},
//                {3374, 2243},
//                {3419, 3401},
//                {2955, 3440},
//                {3150, 3250},
//                {2575, 2257},
//                {2788, 2226},
//                {2310, 2275},
//                {4530, 3212},
//
//                {2304, 2212},
//                {4639, 1315},
//                {2177, 2444},
//                {1712, 1699},
//                {2488, 1735},
//                {5326, 1856},
//                {2238, 1329},
//                {1196, 1404},
//                {2312, 7701},
//                {4386, 5302},
//                {3007, 1470},
//                {2562, 1256},
//                {5788, 1991},
//                {2381, 1076},
//                {4332, 6253},
//                {2715, 1378},
//                {5918, 2579},
//                {2061, 2770},
//                {4780, 2312},
//                {4676, 2478},
//                {1029, 2938},
//                {1263, 2031},
//                {2429, 4308},
//                {2507, 2567},
//                {2394, 3443},
//                {6439, 5501},
//                {4935, 2340},
//                {2140, 4350},
//                {1545, 2557},
//                {1778, 3426},
//                {4370, 5475},
//                {2520, 4312},
//
//                {7374, 2112},
//                {5679, 1215},
//                {6197, 2344},
//                {6772, 1499},
//                {7478, 1635},
//                {8376, 1756},
//                {3278, 1829},
//                {2176, 1204},
//                {1372, 7301},
//                {4376, 5302},
//                {5077, 1570},
//                {6572, 1256},
//                {4778, 1691},
//                {4371, 1576},
//                {5372, 6753},
//                {4775, 1478},
//                {6978, 2579},
//                {5071, 2570},
//                {2770, 2612},
//                {5696, 2378},
//                {6079, 2438},
//                {1273, 2531},
//                {2449, 1708},
//                {4557, 3267},
//                {5374, 2343},
//                {6419, 4301},
//                {7955, 4340},
//                {2150, 5350},
//                {1575, 4457},
//                {4788, 4426},
//                {3310, 3275},
//                {6530, 2312}
//        };

//double city_pos[N][2] =
//{
//                {1304, 2312},
//                {3639, 5315},
//                {4177, 3444},
//                {3712, 5399},
//                {3488, 1335},
//                {3326, 5356},
//                {3238, 3229},
//                {4196, 4304},
//                {4312, 3401},
//                {4386, 5532},
//                {3007, 5370},
//                {2562, 4356},
//                {2788, 3291},
//                {2381, 5376},
//                {1332, 5353},
//                {3715, 6578},
//                {3918, 7579},
//                {4061, 2270},
//                {3780, 2312},
//                {3676, 2578},
//                {4029, 2538},
//                {4263, 5631},
//                {3429, 3408},
//                {3507, 5367},
//                {3394, 6743},
//                {3439, 3901},
//                {2935, 5340},
//                {3140, 2350},
//                {2545, 2357},
//                {2778, 9826},
//                {2370, 9775},
//                {4520, 3912},
//
//                {1374, 1212},
//                {3679, 1315},
//                {4197, 1444},
//                {3772, 1799},
//                {3478, 1435},
//                {3376, 7656},
//                {3278, 1029},
//                {4176, 1904},
//                {4372, 7201},
//                {4376, 1202},
//                {3077, 1430},
//                {2572, 1566},
//                {2778, 1451},
//                {2371, 1326},
//                {1372, 6633},
//                {3775, 1248},
//                {3978, 2539},
//                {4071, 2210},
//                {3770, 2532},
//                {3696, 2688},
//                {4079, 2268},
//                {4273, 2761},
//                {3449, 1968},
//                {3557, 2227},
//                {3374, 2233},
//                {3419, 3351},
//                {2955, 3640},
//                {3150, 3320},
//                {2575, 2527},
//                {2788, 2636},
//                {2310, 2435},
//                {4530, 3532},
//
//                {2304, 2222},
//                {4639, 1655},
//                {2177, 2354},
//                {1712, 1639},
//                {2488, 1765},
//                {5326, 1876},
//                {2238, 1389},
//                {1196, 1494},
//                {2312, 7731},
//                {4386, 5232},
//                {3007, 1210},
//                {2562, 1256},
//                {5788, 1531},
//                {2381, 1536},
//                {4332, 6243},
//                {2715, 1538},
//                {5918, 2539},
//                {2061, 2760},
//                {4780, 2372},
//                {4676, 2868},
//                {1029, 2128},
//                {1263, 2321},
//                {2429, 4348},
//                {2507, 2537},
//                {2394, 3533},
//                {6439, 5231},
//                {4935, 2320},
//                {2140, 4340},
//                {1545, 2567},
//                {1778, 3746},
//                {4370, 5785},
//                {2520, 4362},
//
//                {7374, 2122},
//                {5679, 1245},
//                {6197, 2124},
//                {6772, 1329},
//                {7478, 1435},
//                {8376, 1566},
//                {3278, 1319},
//                {2176, 3204},
//                {1372, 3301},
//                {4376, 5202},
//                {5077, 2470},
//                {6572, 1356},
//                {4778, 6391},
//                {4371, 13576},
//                {5372, 1523},
//                {4775, 6448},
//                {6978, 1279},
//                {5071, 3470},
//                {2770, 2512},
//                {5696, 6478},
//                {6079, 5838},
//                {1273, 7831},
//                {2449, 9708},
//                {4557, 3767},
//                {5374, 4343},
//                {6419, 3601},
//                {7955, 4340},
//                {2150, 5313},
//                {1575, 4454},
//                {4788, 4465},
//                {3310, 3275},
//                {6530, 2311},
//
//                {1754, 2375},
//                {3329, 1335},
//                {4752, 2235},
//                {3322, 1323},
//                {3438, 1553},
//                {3346, 1525},
//                {3658, 1274},
//                {4246, 1045},
//                {4532, 7924},
//                {4736, 5797},
//                {3347, 1954},
//                {2232, 1721},
//                {2238, 1436},
//                {2131, 1634},
//                {1542, 6964},
//                {3345, 1623},
//                {3428, 2153},
//                {4231, 2364},
//                {3640, 2285},
//                {3736, 2586},
//                {4249, 2897},
//                {4533, 2997},
//                {3239, 1965},
//                {3437, 2332},
//                {3324, 2623},
//                {3329, 3234},
//                {2545, 3287},
//                {3340, 3512},
//                {2435, 2323},
//                {2578, 2853},
//                {2640, 2934},
//                {4340, 3464},
//
//                {1094, 2268},
//                {3029, 1287},
//                {4057, 2445},
//                {3072, 1276},
//                {3098, 1254},
//                {3926, 1243},
//                {3948, 1423},
//                {4096, 1276},
//                {4032, 7256},
//                {4576, 5254},
//                {3907, 1287},
//                {2802, 1245},
//                {2708, 1232},
//                {2391, 1254},
//                {1362, 6223},
//                {3795, 1243},
//                {3978, 2254},
//                {4011, 2271},
//                {3450, 2125},
//                {3906, 2237},
//                {4769, 2232},
//                {4273, 2221},
//                {3879, 1257},
//                {3507, 2327},
//                {3724, 2533},
//                {3439, 3251},
//                {1235, 3230},
//                {4320, 3230},
//                {2455, 2257},
//                {2438, 2253},
//                {2120, 2257},
//                {4340, 3234},
//
//                {2254, 2224},
//                {4239, 1353},
//                {2127, 2464},
//                {1742, 1664},
//                {2468, 1775},
//                {5236, 1885},
//                {2348, 1323},
//                {1566, 1432},
//                {2572, 7735},
//                {4436, 5353},
//                {3237, 1475},
//                {2132, 1234},
//                {2338, 1934},
//                {2233, 1073},
//                {4567, 6234},
//                {2432, 1364},
//                {5923, 2586},
//                {2023, 2743},
//                {4753, 2312},
//                {4621, 2423},
//                {1052, 2935},
//                {1232, 2075},
//                {2424, 4312},
//                {2553, 2534},
//                {2364, 3412},
//                {6433, 5535},
//                {4943, 2353},
//                {2163, 4332},
//                {1532, 2543},
//                {1764, 3453},
//                {4323, 5432},
//                {2554, 4387},
//
//                {7374, 8712},
//                {5632, 6715},
//                {6187, 8744},
//                {6754, 9899},
//                {7464, 7635},
//                {8356, 5656},
//                {3232, 5629},
//                {2154, 4704},
//                {1375, 7401},
//                {4332, 6302},
//                {5012, 5470},
//                {6752, 6456},
//                {4723, 4391},
//                {4326, 6476},
//                {5235, 1253},
//                {4752, 3478},
//                {6965, 6479},
//                {5046, 6470},
//                {2764, 4312},
//                {5643, 6478},
//                {6064, 6438},
//                {1264, 6431},
//                {2447, 7408},
//                {4587, 6467},
//                {5375, 4543},
//                {6443, 9801},
//                {7954, 1240},
//                {2155, 2350},
//                {1575, 3457},
//                {4713, 6826},
//                {3332, 3475},
//                {6545, 5612}
//        };

// function statement
double distance(double *, double *);                              // calculate the distance between two cities
double path_len(int *arr, double **city_divided, int CityNumber); // calculate the path length of the TSP solution
void init(int *city_list, int cityNumber);                        // initial function for city list
void create_new(int CityNumber, int *city_list);                  // generate the new solution for TSP

/*
 * distance 
 * DESCRIPTION: calculate the distance between two cities
 * INPUTS:  city1 - first city
 *          city2 - second city
 * OUTPUTS: None
 * RETURN_VALUE: distance between two cities
 * SIDE_EFFECTS: None
 */
double distance(double *city1, double *city2) {

    double x1 = *city1;
    double y1 = *(city1 + 1);
    double x2 = *(city2);
    double y2 = *(city2 + 1);
    double dis = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    return dis;
}

/*
 * path_len 
 * DESCRIPTION: calculate the path length of the TSP solution
 * INPUTS:  arr - the index list for input city in the total city list
 *          city_divided - data coordinates
 *          CityNumber - number of cities in the input city list
 * OUTPUTS: None
 * RETURN_VALUE: The path length on one TSP soluton
 * SIDE_EFFECTS: None
 */
double path_len(int *arr, double **city_divided, int CityNumber) {

    // test
    if (CityNumber == N) {
        double path = 0;  // initial

        // sequential implement
        for (int i = 0; i < CityNumber - 1; i++) {
            int index1 = *(arr + i);
            int index2 = *(arr + i + 1);

            double dis = distance(city_pos[index1 - 1],
                                  city_pos[index2 - 1]);
            path += dis;
        }

        // openMP version
        // #pragma omp parallel
        // {
        //     double temp = 0;
        //     #pragma omp for
        //     for (int i = 0; i < CityNumber - 1; i++)
        //     {
        //         int index1 = *(arr + i);
        //         int index2 = *(arr + i + 1);

        //         double dis = distance(city_pos[index1 - 1],
        //                             city_pos[index2 - 1]);
        //         temp += dis;
        //     }
        //     #pragma omp atomic
        //     path += temp;
        // }

        int last_index = *(arr + CityNumber - 1);   // the index for last city in the global city list
        int first_index = *arr;                     // the index for first city in the global city list
        double last_dis = distance(city_pos[last_index - 1],
                                   city_pos[first_index - 1]);
        path = path + last_dis;
        // return the total length for the solution 
        return path;
    } else {
        double path = 0;


        for (int i = 0; i < CityNumber - 1; i++) {
            int index1 = *(arr + i);
            int index2 = *(arr + i + 1);
            double dis = distance(city_divided[index1 - 1],
                                  city_divided[index2 - 1]);

            path += dis;
        }

        // #pragma omp parallel
        // {
        //     double temp = 0;
        //     #pragma omp for 
        //     for (int i = 0; i < CityNumber - 1; i++)
        //     {
        //         int index1 = *(arr + i);
        //         int index2 = *(arr + i + 1);
        //         double dis = distance(city_divided[index1 - 1],
        //                             city_divided[index2 - 1]);

        //         temp += dis;
        //     }
        //     #pragma omp atomic
        //     path += temp;
        // }

        int last_index = *(arr + CityNumber - 1); // the index for last city in the global city list
        int first_index = *arr;                   // the index for first city in the global city list
        double last_dis = distance(city_divided[last_index - 1],
                                   city_divided[first_index - 1]);
        path = path + last_dis;
        // return total length for the TSP solution
        return path;
    }
}

/*
 * init 
 * DESCRIPTION: initialize one TSP solution for SA to start
 * INPUTS:  City_list - index solution to represent one solution
 *          CityNumber - number of City list
 * OUTPUTS: None
 * RETURN_VALUE: None
 * SIDE_EFFECTS: initialize one TSP solution for SA to start
 */

void init(int *city_list, int cityNumber) {
    for (int i = 0; i < cityNumber; i++)
        city_list[i] = i + 1; // generate the initial solution for TSP
}

/*
 * create_new 
 * DESCRIPTION: create a new TSP solution for TS algro.
 * INPUTS:  City_list - index solution to represent one solution
 *          CityNumber - number of City list
 * OUTPUTS: None
 * RETURN_VALUE: None
 * SIDE_EFFECTS: create a new TSP solution for TS algro. by randomly swapping two points.
 */

void create_new(int CityNumber, int *city_list) {
    double r1 = ((double) rand()) / (RAND_MAX + 1.0);
    double r2 = ((double) rand()) / (RAND_MAX + 1.0);
    int pos1 = (int) (CityNumber * r1);
    int pos2 = (int) (CityNumber * r2);

    // swap two points
    int temp = city_list[pos1];
    city_list[pos1] = city_list[pos2];
    city_list[pos2] = temp;
}

/*
 * swap 
 * DESCRIPTION: swap two elements in the array, a[i] and a[offset]
 * INPUTS:  a - array
 *          i - one position 
 *          offset - second position
 * OUTPUTS: None
 * RETURN_VALUE: None
 * SIDE_EFFECTS: swap two elements in the array, a[i] and a[offset]
 */
void swap(int i, int *a, int offset) {
    int temp;
    temp = a[offset];
    a[offset] = a[i];
    a[i] = temp;
}

/*
 * perm 
 * DESCRIPTION: return the global TSP solution and its cities id
 * INPUT: midlist - result index sequence for global
 *        citylist - receive list from other processes
 *        index - the subscript for city lists from different processes
 *        len - number of iteration
 *        remain - remain cities
 *        offset  - offset for swapping index
 * OUTPUTS: midlist - the id of cities for the solution
 * RETURN_VALUE: shortest path length
 * SIDE_EFFECTS: using the method of Full arrangement to generate the TSP solution from receving partial solution of other processes
 *               and return shortest path length
 */
double perm(int *midlist, int **citylist, int *index, int len, int cityNumber, int remain, int offset) {
    // len is the iteration number
    // nn records the first index of midlist
    int n = 0; // record the second index of midlist
    static int nn = 0;
    static double shortestListDistance;  // return value
    double distance;
    // int list[N];
    int *list = (int *) malloc(N * sizeof(int));
    if (offset == len - 1) {
        if (remain == 0) {
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < cityNumber; j++) {
                    list[n] = citylist[index[i]][j];
                    // printf("nn=%d  , n=%d  List=%d\n",nn,n,list[n]);
                    n++;
                }
            }
        } else {
            for (int i = 0; i < len; i++) {
                if (index[i] <= (remain - 1)) {
                    for (int j = 0; j < cityNumber + 1; j++) {
                        list[n] = citylist[index[i]][j];
                        n++;
                    }
                } else {
                    for (int j = 0; j < cityNumber; j++) {
                        // printf("n=%d",n);
                        list[n] = citylist[index[i]][j];
                        n++;
                    }
                }
            }
        }

        // check whether it is the shortest list distance
        if (nn == 0) {
            shortestListDistance = path_len(list, (double **) city_pos, N);
            for (int i = 0; i < N; i++) {
                midlist[i] = list[i];
            }
        } else {
            distance = path_len(list, (double **) city_pos, N);
            if (distance <= shortestListDistance) {
                shortestListDistance = distance;
                for (int i = 0; i < N; i++) {
                    midlist[i] = list[i];
                }
            }
        }
        nn++;
        return shortestListDistance;
    }
    for (int i = offset; i < len; i++) {
        swap(i, index, offset);
        free(list);
        perm(midlist, citylist, index, len, cityNumber, remain, offset + 1);
        swap(i, index, offset);
    }
    return shortestListDistance;
}

// main 
int main(int argc, char *argv[]) {
    // cityNumber is the number of cities each process needs to handle with
    int myid, numprocs, namelen, cityNumber;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    double **city_divided;
    double startwtime = 0.0, midwtime, endwtime;

    // MPI set
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);


    if (numprocs > 1) {

        // divide 
        int remain = N % (numprocs - 1);
        if (remain != 0) {
            if (myid <= remain && myid != 0) {

                cityNumber = ((N - remain) / (numprocs - 1)) + 1;
                city_divided = (double **) malloc(cityNumber * sizeof(double));
                for (int i = 0; i < cityNumber; i++) {
                    city_divided[i] = (double *) malloc(2 * sizeof(double));
                }
                for (int i = 0; i < cityNumber; i++) {
                    // printf("test,myid=%d\n",myid);
                    city_divided[i][0] = city_pos[cityNumber * (myid - 1) + i][0];
                    city_divided[i][1] = city_pos[cityNumber * (myid - 1) + i][1];
                }
            } else {
                cityNumber = ((N - remain) / (numprocs - 1));
                city_divided = (double **) malloc(cityNumber * sizeof(double));
                for (int i = 0; i < cityNumber; i++) {
                    city_divided[i] = (double *) malloc(2 * sizeof(double));
                }
                for (int i = 0; i < cityNumber; i++) {
                    city_divided[i][0] = city_pos[cityNumber * (myid - 1) + remain + i][0];
                    city_divided[i][1] = city_pos[cityNumber * (myid - 1) + remain + i][1];
                }
            }
        } else {
            cityNumber = N / (numprocs - 1);
            city_divided = (double **) malloc(cityNumber * sizeof(double));
            for (int i = 0; i < cityNumber; i++) {
                city_divided[i] = (double *) malloc(2 * sizeof(double));
            }
            for (int i = 0; i < cityNumber; i++) {
                city_divided[i][0] = city_pos[cityNumber * (myid - 1) + i][0];
                city_divided[i][1] = city_pos[cityNumber * (myid - 1) + i][1];
            }
        }

        if (myid == 0) {
            // cityList is used to store the result of sub processes
            int **cityList = (int **) malloc((numprocs - 1) * sizeof(int));

            if (remain == 0) {
                for (int i = 0; i < numprocs - 1; i++) {
                    cityList[i] = (int *) malloc(cityNumber * sizeof(int));
                }
            } else {
                for (int i = 0; i < numprocs - 1; i++) {

                    if (i < remain) {
                        cityList[i] = (int *) malloc((cityNumber + 1) * sizeof(int));
                    } else {
                        cityList[i] = (int *) malloc(cityNumber * sizeof(int));
                    }
                }
            }

            startwtime = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);
            midwtime = MPI_Wtime();
            // cTime keeps the calculation time of sub nodes
            double cTime = midwtime - startwtime;
            printf("slaver finished,time:%f s \t ,final linking...\n", cTime);
            // start to receive all results
            if (remain == 0) {
                for (int i = 0; i < numprocs - 1; i++) {

                    MPI_Recv(cityList[i], cityNumber, MPI_INT, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } else {
                for (int i = 0; i < numprocs - 1; i++) {
                    if (i <= remain) {
                        MPI_Recv(cityList[i], cityNumber + 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } else {
                        MPI_Recv(cityList[i], cityNumber, MPI_INT, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            // start to  connect routines
            int index[numprocs - 1];
            for (int i = 0; i < numprocs - 1; i++) {
                index[i] = i;
            }

            int allpai = numprocs - 1;
            for (int i = numprocs - 2; i > 0; i--) {
                allpai = allpai * i;
            }

            // bestList keeps local optimal routine
            int bestList[N];
            double shortestListDistance1;
            shortestListDistance1 = perm(bestList, cityList, index, numprocs - 1, cityNumber, remain, 0);

            endwtime = MPI_Wtime();

            printf("\n\nall finished,shortest distance is :%f\t,total time :%f s\t\n", shortestListDistance1,
                   (endwtime - startwtime));
            printf("relative best path :\n");

            // output the optimal routine
            for (int i = 0; i < N; i++) {
                if (i == N - 1) {
                    printf("%d \n", bestList[i]);
                } else {
                    printf("%d--->", bestList[i]);
                }
            }
        } else {

            int city_list[cityNumber];

            double T;

            // record annealing time 
            int count = 0;
            // initial temperature 
            T = T0;
            // initialize a solution
            init(city_list, cityNumber);

            // use city_list_copy to save the original solution
            int city_list_copy[cityNumber];

            // f1 is the inital solution's corresonding object function value
            double f1;
            // f2 is the generated solution's corresponding object function value
            double f2;
            // df is the difference betwen f1 and f2
            double df;


            // a generated random number between 0 and 1
            double r;

            // stop simulated annealing if the temperature is lower than the target temperature
            while (T > T_end) {
                for (int i = 0; i < L; i++) {
                    // copy list
                    memcpy(city_list_copy, city_list, cityNumber * sizeof(int));
                    // generate a new solution
                    create_new(cityNumber, city_list);
                    f1 = path_len(city_list_copy, city_divided, cityNumber);
                    f2 = path_len(city_list, city_divided, cityNumber);

                    df = f2 - f1;
                    // apply Metropolis principle
                    if (df >= 0) {
                        r = ((double) rand()) / (RAND_MAX);

                        // if the newly generated solution is worse, reverse back to initial solution
                        if (exp(-df / T) <= r) {
                            memcpy(city_list, city_list_copy, cityNumber * sizeof(int));
                        }
                    }
                }
                // Anneal
                T *= q;
                count++;
            }

            // transfer local index into global index 
            for (int i = 0; i < cityNumber; i++) {
                if (myid <= remain) {
                    city_list[i] += ((myid - 1) * cityNumber);
                } else {
                    city_list[i] += ((myid - 1) * cityNumber + remain);
                }
            }

            printf("No.%d proc which is on the mechine %s finished , the path is:\n", myid, processor_name);
            for (int i = 0; i < cityNumber; i++) {
                if (i == cityNumber - 1) {
                    printf("%d \n", city_list[i]);
                } else {
                    printf("%d -->", city_list[i]);
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Send(city_list, cityNumber, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    } else // case for one processor
    {
        printf(" run the sernial version of tsp \n");
        int city_list[N];
        // record annealing time 
        int count = 0;
        // initial temperature 
        double T = T0;
        // // initialize a solution
        init(city_list, N);
        // use city_list_copy to save the original solution
        int city_list_copy[N];
        // f1 is the inital solution's corresonding object function value
        double f1;
        // f2 is the generated solution's corresponding object function value
        double f2;
        // df is the difference betwen f1 and f2
        double df;

        // random number between 0 and 1
        double r;

        startwtime = MPI_Wtime();

        // stop simulated annealing if the temperature is lower than the target temperature
        while (T > T_end) {
            for (int i = 0; i < L; i++) {
                memcpy(city_list_copy, city_list, N * sizeof(int));
                // generate new solution
                create_new(N, city_list);
                f1 = path_len(city_list_copy, (double **) city_pos, N);
                f2 = path_len(city_list, (double **) city_pos, N);
                df = f2 - f1;
                // Metropolis principle
                if (df >= 0) {
                    r = ((double) rand()) / (RAND_MAX);
                    // if the newly generated solution is worse, reverse back to initial solution
                    if (exp(-df / T) <= r) {
                        memcpy(city_list, city_list_copy, N * sizeof(int));
                    }
                }
            }
            T *= q; // annealing
            count++;
        }

        endwtime = MPI_Wtime();

        float shortestListDistance1 = path_len(city_list_copy, (double **) city_pos, N);
        printf("\n\nall finished,shortest distance is :%f\t,total time :%f s\t\n", shortestListDistance1,
               (endwtime - startwtime));
        printf("relative best path :\n");
        // print out the optimal routine
        for (int i = 0; i < N; i++) {
            if (i == N - 1) {
                printf("%d \n", city_list[i]);
            } else {
                printf("%d--->", city_list[i]);
            }
        }
    }

    MPI_Finalize();

    return 0;
}

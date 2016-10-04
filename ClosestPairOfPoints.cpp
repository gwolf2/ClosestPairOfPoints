/********************************************************************************************
Name: Gavin Wolf		Z#: 15289719
Course: COT 6405: Analysis of Algorithms
Professor: Dr. Mihaela Cardei
Due Date: 4/21/2016
Programming Project

Description: This program implements two algorithms to compute the closest pair of points
 in a set, and measures their actual run-times over a number of trials.
*********************************************************************************************/

#include <iostream>
#include <float.h>
#include <random>
using namespace std;

/********************************************************************************************
User inputs.
*********************************************************************************************/

//Enter input size up to 20000 to run:
int g = 20000;

//Select algorithm:
// 1 for bruteForce
// 2 for closestPair (divide and conquer)
int select = 2;

//Enter number of trials
int trials = 20;

/********************************************************************************************
Struct, function and array declarations.
*********************************************************************************************/

//A structure to represent a point in a plane
struct Point
{
    unsigned int x, y;
};

//Function declarations
double dist(Point p1, Point p2);
double bruteForce(Point P[], int n);
double SClosest(Point *S, int size, double d);
double closestPairRec(Point *Px, Point *Py, int n);
double closestPair(Point *P, int n);
void Merge(Point *M, int p, int q, int r, int compareOption);
void MergeSort(Point *M, int p, int r, int compareOption);

//Array declarations
Point P[20000];

// for MergeSort
Point L[10001]; //+1 to account for sentinel value
Point R[10001]; //+1 to account for sentinel value

/********************************************************************************************
Simulation and running time measurements.
*********************************************************************************************/
int main()
{
    //Initialize seed for random number generator
    srand(time(NULL));

    //Heading for output
    cout << endl << "Results:" << endl << endl << "# Elements\t\tRun Time\t\tTrials" << endl;

    int sumTime = 0; //Will be used to calculate the average

    //Loop to run multiple trials
    for (int i = 0; i < trials; i++)
    {
        //Fill array P with random values.
        // Note: rand() + rand() used instead of just rand() to increase the likelihood that
        //  the arrays will have distinct values.
        // Note: It is important to generate a new set of random numbers for each trial in
        //  order to test the algorithms' run times on different inputs.
        for (int j = 0; j < g; j++)
        {
            P[j].x = (unsigned int)(rand() + rand() - 1);
            P[j].y = (unsigned int)(rand() + rand() - 1);
        }

        //Store startTime of algorithm execution
        auto startTime = std::chrono::high_resolution_clock::now();

        //Execute algorithm
        if (select == 1)
        {
            bruteForce(P, g);
        }
        else
        {
            closestPair(P, g);
        }

        //Store finishTime of algorithm execution and calculate runTime
        auto finishTime = std::chrono::high_resolution_clock::now();
        auto runTime = finishTime - startTime; //Execution time
        long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(runTime).count();

        sumTime += microseconds; //When loop exits, this will be the sum of all the trials
    }

    //Calculate average run time over all trials
    int averageTime = sumTime / trials;

    //Output results
    cout << g << "\t\t\t" << averageTime << "\t\t\t\t" << trials << endl;

    return 0;
}

/********************************************************************************************
closestPair function - Divide and conquer algorithm
*********************************************************************************************/
//Constructs Px (all points sorted by x-coordinate in increasing order)
//Constructs Py (all points sorted by y-coordinate in increasing order)
//Calls closestPairRec (recursive function) on Px/Py
//Returns the closest pair of points
double closestPair(Point *P, int n)
{
    Point Px[n];
    Point Py[n];
    for (int i = 0; i < n; i++)
    {
        Px[i] = P[i];
        Py[i] = P[i];
    }

    //Use MergeSort to sort arrays
    MergeSort(Px, 0, n - 1, 1); //1 for compare x
    MergeSort(Py, 0, n - 1, 2); //2 for compare y

    //Call closestPairRec to find the smallest distance
    return closestPairRec(Px, Py, n);
}

/********************************************************************************************
closestPairRec function - Recursive function
*********************************************************************************************/
//Recursive function to find the smallest distance
double closestPairRec(Point *Px, Point *Py, int n)
{
    //If there are 2 or 3 points, just use brute force
    if (n <= 3)
    {
        return bruteForce(Px, n);
    }

    //Find the middle point of the x-coordinates
    int mid = ((n - 1) / 2);
    Point midPoint = Px[mid];

    //Divide points in y sorted array around the vertical line.
    Point Qy[mid + 1]; //The left-half of P, sorted by y-coordinate
    Point Ry[n - mid - 1]; //The right-half of P, sorted by y-coordinate

    //Indexes of left and right subarrays
    int leftIndex = 0;
    int rightIndex = 0;

    //Construct Qy and Ry by traversing through Py
    for (int i = 0; i < n; i++)
    {
        if (Py[i].x <= midPoint.x)
        {
            Qy[leftIndex] = Py[i];
            leftIndex++;
        }
        else
        {
            Ry[rightIndex] = Py[i];
            rightIndex++;
        }
    }

    //Note: the pseudocode in Kleinberg & Tardos indicates to construct the following arrays:
    // Qx - The left-half of P, sorted by x-coordinate
    // Rx - The right-half of P, sorted by x-coordinate
    //Instead of constructing these arrays explicitly, the code below uses the arrays
    // implicitly by using "mid" to delimit the left and right halves of P.

    //Recursively find the closest pair among the points in Q and R
    double dl = closestPairRec(Px, Qy, mid);
    double dr = closestPairRec(Px + mid, Ry, n - mid - 1);

    //Take the minimum of the two distances
    double d = (dl < dr) ? dl : dr;

    //S is is the set of points in P within d distance of a line passing through the
    // midpoint of the x-coordinates, sorted by y-coordinate.
    //Note: S has the property that if any two points within it are closer to each other
    // than d, then those two points are within 15 positions of each other in S.
    // Therefore, the closest pair of points in S can be computed in linear time,
    // faster than the quadratic time found in the brute force method of finding the
    // closest pair of points in a set.
    Point S[n];
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        long long int abso = 0;

        if ((Py[i].x - midPoint.x) < 0)
        {
            abso = -(Py[i].x - midPoint.x);
        }
        else
        {
            abso = (Py[i].x - midPoint.x);
        }

        if (abso < d)
        {
            S[j] = Py[i];
            j++;
        }
    }

    //Call helper function SClosest to computer the minimum distance in S
    return SClosest(S, j, d);
}

/********************************************************************************************
SClosest function
*********************************************************************************************/
//Helper function to find the distance between the closest pair of points in S.
double SClosest(Point *S, int size, double d)
{
    double min = d; // Initialize the minimum distance as d

    //If two points are found whose y coordinates are closer than the minimum distance,
    // then investigate whether those points are closer than the minimum distance, and keep
    // track of the minimum value so far.
    for (int i = 0; i < size; i++)
    {
        for (int j = i + 1; j < size && (S[j].y - S[i].y) < min; j++)
        {
            if (dist(S[i], S[j]) < min)
            {
                min = dist(S[i], S[j]);
            }
        }
    }

    return min;
}

/********************************************************************************************
bruteForce function
*********************************************************************************************/
//Computes the closest pair of points in P in a brute force manner - by computing the
// distance each pair of points and returning the minimum distance.
double bruteForce(Point *P, int n)
{
    double min = DBL_MAX;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (dist(P[i], P[j]) < min)
            {
                min = dist(P[i], P[j]);
            }
        }
    }
    return min;
}

/********************************************************************************************
dist function
*********************************************************************************************/
//Computes the Euclidean distance between two points.
double dist(Point p1, Point p2)
{
    //Casting each int value to long long int so that the values can be safely multiplied,
    // i.e., without the risk of "overflowing" the capacity of int.
    long long int p1x = p1.x;
    long long int p2x = p2.x;
    long long int p1y = p1.y;
    long long int p2y = p2.y;

    return sqrt( (p1x - p2x) * (p1x - p2x) + (p1y - p2y) * (p1y - p2y) );
}

/********************************************************************************************
MergeSort function
*********************************************************************************************/
//Recursively sorts the array of structs M by either the x-coordinates, or the y-coordinates.
// If compareOption == 1, the array of structs will be sorted based on the x-coordinate;
// otherwise, the array of structs will be sorted based on the y-coordinate.
void MergeSort(Point *M, int p, int r, int compareOption)
{
    if (p < r)
    {
        int q = (p + r) / 2;
        MergeSort(M, p, q, compareOption);
        MergeSort(M, q + 1, r, compareOption);
        Merge(M, p, q, r, compareOption);
    }
}

/********************************************************************************************
Merge function
*********************************************************************************************/
//A helper function for MergeSort
void Merge(Point *M, int p, int q, int r, int compareOption)
{
    int n1 = q - p; //Last index of left array
    int n2 = r - q - 1; //Last index of right array

    //Construct the left half of array M
    for (int i = 0; i <= n1; i++)
        L[i] = M[p + i];

    //Construct the right half of array M
    for (int j = 0; j <= n2; j++)
        R[j] = M[q + j + 1];

    //Depending on the value of compareOption, combine the left and right halves of array
    // back into M.
    if (compareOption == 1)
    {
        L[n1 + 1].x = RAND_MAX + RAND_MAX; //The sentinel value
        R[n2 + 1].x = RAND_MAX + RAND_MAX; //The sentinel value

        int i = 0;
        int j = 0;

        for (int k = p; k <= r; k++)
        {
            if (L[i].x <= R[j].x)
            {
                M[k] = L[i];
                i++;
            }
            else
            {
                M[k] = R[j];
                j++;
            }
        }
    }
    else { //sort based on Y-coordinate
        L[n1 + 1].y = RAND_MAX + RAND_MAX; //The sentinel value
        R[n2 + 1].y = RAND_MAX + RAND_MAX; //The sentinel value

        int i = 0;
        int j = 0;

        for (int k = p; k <= r; k++)
        {
            if (L[i].y <= R[j].y)
            {
                M[k] = L[i];
                i++;
            }
            else
            {
                M[k] = R[j];
                j++;
            }
        }
    }
}

/*
***SAMPLE OUTPUT***

Results:

# Elements		Run Time		Trials
20000			14828				20

*/
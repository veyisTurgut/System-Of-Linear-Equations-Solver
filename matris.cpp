#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h> /* pow */

void reader(std::ifstream &stream);
void rrefToIdentity(int n, double *matris, std::string type);
void findArbitraryVars(int n, double *matris);
void solver(int n, double *matris);
void swapRows(int n, double *row1, double *row2);
void multiplyByScalar(int n, double *row, double scalar);
void subtractRowWithAMultiplier(int n, double *row1, double *row2, double scalar);
int main(int argc, char *argv[])
{

    std::ifstream inFile;
    for (int i = 1; i < argc; i++) //read the arguments
    {
        inFile.open(argv[i]);
        if (!inFile)
        {
            std::cout << "Unable to open file " << argv[i];
            exit(1); // terminate with error
        }
        std::cout << "-------------------------------------------\n";
        reader(inFile);
        inFile.close();
    }
    std::cout << "-------------------------------------------\n";
    return 0;
}

void reader(std::ifstream &stream)
{
    int n;
    stream >> n;
    double matris[n][2 * n + 1];
    for (int i = 0; i < n; i++) //read the file
    {
        for (int j = 0; j < 2 * n + 1; j++)
        {
            if (j < n + 1)
                stream >> matris[i][j];
            else if (j - i == n + 1)
                matris[i][j] = 1;
            else
                matris[i][j] = 0;
        }
    }
    //now matris = A|b|I
    solver(n, (double *)matris);
}

void solver(int n, double *matris)
{
    //*((arr+i*n) + j) ile arr[i][j] ye eriş.
    double pivot;
    int indexOfPivot;
    for (int i = 0; i < n; i++) //start the algorithm //n represents number of rows //i represents current row
    {

        indexOfPivot = i * (2 * n + 2);
        pivot = *(matris + indexOfPivot); // choose pivot
        if (abs(pivot) < 0.00001)         //pivot = 0, choose another pivot
        {
            for (int j = i + 1; j < n; j++) //find next nonzero element in this column.
            {
                double candidatePivot = *(matris + j * (2 * n + 1) + i); //searching elements under current pivot
                if (abs(candidatePivot) > 0.00000001)                    //if its nonzero, swap rows
                {
                    swapRows(n, (double *)matris + i * (2 * n + 1), (double *)matris + j * (2 * n + 1));
                    break; //already found the correct pivot, no need to continue.
                }
            } //now we swapped rows and pivot is nonzero( hopefully)
        }
        pivot = *(matris + indexOfPivot);
        if (abs(pivot) < 0.000001)
        { //means we couldn't swap.
            if (abs(*(matris + n + i * (2 * n + 1))) < 0.0000001)
                findArbitraryVars(n, (double *)matris);
            else
                std::cout << "\nInconsistent problem\n\n";
            return;
        }
        //normalize the row which contains pivot.
        multiplyByScalar(n, (double *)matris + i * (2 * n + 1), 1 / pivot);

        for (int j = i + 1; j < n; j++)
        { //multiply each rows with (1/elemnt under pivot)

            multiplyByScalar(n, (double *)matris + j * (2 * n + 1), 1 / (*(matris + indexOfPivot + (j - i) * (2 * n + 1))));
            //normalized rows
            //now it's time to subtract other rows from pivot row so that only pivot is 1 in its column
            if (abs(*(matris + indexOfPivot + (j - i) * (2 * n + 1))) > 0.000001)
                subtractRowWithAMultiplier(n, (double *)matris + i * (2 * n + 1), (double *)matris + j * (2 * n + 1), 1.0);
        }
    }
    std::cout << "\n";
    rrefToIdentity(n, (double *)matris, "Unique");
    std::cout << "\n";
}

/**
 * Interchange given rows
 */
void swapRows(int n, double *row1, double *row2)
{
    double temp;
    for (int i = 0; i < 1 + 2 * n; i++)
    {
        temp = *(row1 + i);
        *(row1 + i) = *(row2 + i);
        *(row2 + i) = temp;
    }
}

/**
 * Multiply each element in given row by given scalar
 */
void multiplyByScalar(int n, double *row, double scalar)
{
    if (abs(1 / scalar) < 0.000001)
        return;
    for (int i = 0; i < 1 + 2 * n; i++)
    {
        *(row + i) *= scalar;
    }
}

/**
 * row2 = row2 - row1 * scalar
 */
void subtractRowWithAMultiplier(int n, double *row1, double *row2, double scalar)
{
    for (int i = 0; i < 2 * n + 1; i++)
    {
        *(row2 + i) -= *(row1 + i) * scalar;
    }
}
/** 
* If all pivots are present, then this method handles the rest.
* it converts Row reduced echelon form of initial matrix 
* to Identity matrix so that ıdentity matrix becomes A'-1
* If type is "Unique" then prints inverted matrix, else type is 
* "arbitrary" then prints only solution.
*/
void rrefToIdentity(int n, double *matris, std::string type)
{
    for (int i = n - 2; i >= 0; i--)
    {
        for (int j = i; j >= 0; j--)
        {
            //row j -= row i+1 * m[j][i+1]
            subtractRowWithAMultiplier(n, (double *)matris + (i + 1) * (2 * n + 1), (double *)matris + j * (2 * n + 1), *(matris + j * (2 * n + 1) + i + 1));
        }
    }
    if (!type.compare("Unique")) 
        std::cout << "Unique solution: ";
    for (int i = 0; i < n; i++)
    {
        std::cout << *(matris + i * (2 * n + 1) + n) << "  ";
    }
    if (!type.compare("Unique"))
    {
        std::cout << "\nInverted A :\t";
        for (int i = 0; i < n; i++)
        {
            if (i > 0)
                std::cout << "\t\t";
            for (int j = 0; j < n; j++)
            {
                std::cout <<std::setprecision(5)<< *(matris + i * (2 * n + 1) + j + n + 1) << "\t|";
            }
            std::cout << "\n";
        }
    }
}
/**
 * When one (or more) of the pivots is (are) not present, this method handles the rest. 
 * Gives arbitrary variable 0 as value and call unique solution function.
 */ 
void findArbitraryVars(int n, double *matris)
{
    std::cout << "\nArbitrary variables: "; //<<var;
    for (int i = 0; i < n; i++)
    {
        if (abs(*(matris + i * (2 * n + 1) + i) < 0.000001))
        {
            std::cout << "x" << i + 1 << " ";
            //make all elements ,which shares the same colums with this arbitrary variable, 0
            for (int j = 0; j < n; j++)
            {
                *(matris + j * (2 * n + 1) + i) = 0;
            }
        }
    }

    std::cout << "\nArbitrary solution: ";
    rrefToIdentity(n, (double *)matris, "Arbitrary");
    std::cout << "\n\n";
}
#include <iostream>
#include <vector>
using namespace std;

void multiply_matrices(vector <vector<int> > matrix1,vector <vector<int> > matrix2, int cols, int rows2); 
void setMatrix(vector <vector<int> > &matrix, int rows, int cols);

int main()
{
    int rows, cols, rows2, cols2;    
    vector< vector<int> > matrix, matrix2;        
    cout<<"Please enter the number of Rows and Columns for your first Matrix."<<endl;
    cout<<"Rows: ";
    cin>>rows;
    cout<<"Columns: ";
    cin>>cols;

    matrix.resize(cols, vector<int>(rows,0));  //Saw this online so not sure how it works but it works, if i take out one i cant do row<column and vice versa
    matrix.resize(rows, vector<int>(cols,0));

    cout<<"Size has been declared, please enter data for your matrix"<<endl;

    setMatrix(matrix,rows,cols);

    cout<<"Second Matrix Automatically Set by Matrix Multiplication Rule"<<endl; //Just automatically sets second matrix as per Matrix Multiplication Rule
    rows2=cols;
    cols2=rows;
    cout<<"Second Matrix Size is: " << rows2 << " by " << cols2 << endl;
    matrix2.resize(cols2, vector<int>(rows2,0));
    matrix2.resize(rows2, vector<int>(cols2,0));

    setMatrix(matrix2,rows2,cols2);        

    cout<<"Multiplied Matrix is:"<<endl;
    multiply_matrices(matrix,matrix2,cols,rows2);

    system("PAUSE");
    return 0;
}

void setMatrix(vector <vector<int> > &matrix, int rows,int cols){
     int num;
     for(int i = 0; i < rows; i ++)
     {
        for (int j = 0; j < cols; j++)
        {
            cout << "Enter Value for Row " << (i+1) << " Column " << (j+1) << ": ";
            cin>>num;
            matrix[i][j]=num;            
        }        
        cout <<  endl;
    }

 /*for(int i = 0; i < rows; i ++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << matrix[i][j] << " ";
        }       
        cout <<  endl;
    }          
  */   
     }
void multiply_matrices(vector <vector<int> > matrix1,vector <vector<int> > matrix2, int cols, int rows2){
    vector< vector<int> > tempMatrix;
    int newrows=rows2;
    int newcols=cols;
    int sum;
    tempMatrix.resize(newcols, vector<int>(newrows,0));   //Resizing new matrix to proper size, so if it was (2x3)(3x2), new matrix is (3x3)

    for (int i = 0; i < newrows; i++)                    //This Works Fine for Square Matrixes but not for others, i have no clue how to fix it?
    {
        for (int j = 0; j < newcols; j++){
            //sum=0;    
            for (int u = 0; u < newcols; u++)
            {
                //sum+=matrix1[i][u] * matrix2[u][j];
                //tempMatrix[i][j]=sum;
                tempMatrix[i][j] += matrix1[i][u] * matrix2[u][j];
            }
        }
    }
    for(int i = 0; i < newrows; i ++)
    {
        for (int j = 0; j < newcols; j++)
        {
            cout << tempMatrix[i][j] << " ";
        }        
        cout <<  endl;
    }          
}
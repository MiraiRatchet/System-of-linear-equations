#ifndef CONSTANTS_H
#define CONSTANTS_H
const double epsilon = 1e-4;
const int maxLength = 64;
#else
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;
typedef double T;

bool checkConvergenceJacobiSeidel(const int size, T** matrix) {
    T sum;
    for (int i = 0; i < size; i++) {
        sum = 0;
        for (int j = 0; j < i-1; j++) {
            sum += abs(matrix[i][j]);
        }
        for (int j = i+1; j < size; j++) {
            sum += abs(matrix[i][j]);
        }
        if (sum >= abs(matrix[i][i])) {
            return false;
        }
    }
    return true;
}

T vectorNorm(T* vector, const int size) {
    T norm = 0;
    for (int i = 0; i < size; i++) {
        norm += vector[i] * vector[i];
    }
    return sqrt(norm);
}

T matrixNorm(const int size, T** matrix, T& A1, T& Ainf) {
    T* vertical = new T[size];
    T* horizontal = new T[size];

    for (int i = 0; i < size; i++) {
        vertical[i]=0;
        horizontal[i]=0;
        for (int j = 0; j < size; j++) {
            vertical[i] += abs(matrix[j][i]);
            horizontal[i] += abs(matrix[i][j]);
        }
    }

    A1 = vertical[0];
    Ainf = horizontal[0];

    for (int i = 1; i < size; i++) {
        if (vertical[i] > A1) {
            A1 = vertical[i];
        }
        if (horizontal[i] > Ainf) {
            Ainf = horizontal[i];
        }
    }
    delete[]vertical;
    delete[]horizontal;
    return Ainf;
    //return A1;
}

void matrixCNormJZ(const int size, T** matrix, T** C, T& normC) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i != j) {
                C[i][j] = -matrix[i][j] / matrix[i][i];
            }
            else {
                C[i][j] = 0;
            }
        }
    }
    T Cinf;
    matrixNorm(size, C, normC, Cinf);
    cout<<"Norm 1 of matrix C in Jacobi and Seidel methods: ";
    cout<<normC<<endl;
    cout<<"Norm inf of matrix C in Jacobi and Seidel methods: ";
    cout<<Cinf<<endl;
}

bool checkSimpleIterConvergence(const int size, T** matr, T** C,T tau,T& normC) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                C[i][i] = -(tau*matr[i][i] - 1);
            }
            else {
                C[i][j] = -matr[i][j];
            }
        }
    }
    T Cinf;
    matrixNorm(size, C, normC, Cinf);
    cout<<"Norm 1 of matrix C in simple iterations method: ";
    cout<<normC<<endl;
    cout<<"Norm inf of matrix C in simple iterations method: ";
    cout<<Cinf<<endl;
    return 0;
}

void outputMatrix(T** matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

T residual(T** matrix, T* vector, T* answer, const int size) {
    T* multiplication = new T[size];
    T* difference = new T[size];
    T norm = 0;
    for (int i = 0; i < size; i++) {
        multiplication[i] = 0;
        for (int j = 0; j < size; j++) {
            multiplication[i] += matrix[i][j] * vector[j];
        }
    }
    for (int i = 0; i < size; i++) {
        difference[i] = answer[i] - multiplication[i];
    }
    for (int i = 0; i < size; i++) {
        norm += difference[i]*difference[i];
    }
    norm = sqrt(norm);
    delete[]multiplication;
    delete[]difference;
    return norm;
}

void multiplyMatrices(T** matr1, T** matr2, T** matr3, const int size) {
    for (int a = 0; a < size; a++) {
        for (int b = 0; b < size; b++) {
            matr1[a][b] = 0;
            for (int k = 0; k < size; k++) {
                matr1[a][b] += matr2[a][k] * matr3[k][b];
            }
        }
    }
}

int readSettings(const char* const fileName, size_t& size, char* const matrixName, char* const vectorName) {
    ifstream ifile;
    ifile.open(fileName);
    if (!ifile.is_open()) {
        cerr << " Error: file with settings is not open!\n";
        return 1;
    }
    ifile >> size;
    string str;
    getline(ifile, str);
    ifile >> matrixName;
    getline(ifile, str);
    ifile >> vectorName;
    ifile.close();
    return 0;
}

int readFileMatrixAndVector(char* const fileName1, char* const fileName2, T** matrix, T* vector, const int size) {
    ifstream ifile;
    ifile.open(fileName1);
    if (!ifile.is_open()) {
        cerr << " Error: file cannot be opened!\n";
        return 1;
    }

    cout << "Matrix A: " << endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            ifile >> matrix[i][j];
            cout << matrix[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
    ifile.close();

    ifile.open(fileName2);
    if (!ifile.is_open()) {
        cerr << " Error: file cannot be opened!\n";
        return 1;
    }

    cout << "Vector b: " << endl;
    for (int i = 0; i < size; i++) {
        ifile >> vector[i];
        cout << vector[i] << endl;
    }
    cout << endl;
    ifile.close();

    return 0;
}

int JacobiMethod(const int size, T** matrix, T* vector) {
    T normC = 0;
    T* answer = new T[size];
    T* temp = new T[size];
    for (int i = 0; i < size; i++) {
        answer[i] = 0;
        temp[i] = 0;
    }
    T sum1;
    T sum2;
    T iter = 0;
    T A1 = 0;
    T Ainf = 0;
    T** C = new T* [size];
    for (int i = 0; i < size; i++) {
        C[i] = new T[size];
    }
    T* difference = new T[size];
    for (int i = 0; i < size; i++) {
        difference[i] = 0;
    }
    matrixCNormJZ(size, matrix, C, normC);
    do {
        iter++;
        for (int i = 0; i < size; i++) {
            sum1 = 0;
            sum2 = 0;
            for (int j = i+1; j < size; j++) {
                sum2 += matrix[i][j] * temp[j];
            }
                for (int l = 0; l < i; l++) {
                    sum1 += matrix[i][l] * temp[l];
                }
            answer[i] = (-sum1 - sum2 + vector[i])/ matrix[i][i];
        }
        for (int k = 0; k < size; k++) {
            difference[k] = answer[k] - temp[k];
        }
        for (int k = 0; k < size; k++) {
            temp[k] = answer[k];
        }
    } while ((residual(matrix,answer,vector,size)/matrixNorm(size,matrix,A1,Ainf))>=epsilon);
    cout << "Solved by Jacobi method: " << endl;
    for (int i = 0; i < size; i++) {
        cout << answer[i] << endl;
    }
    cout<<iter<<" iterations for Jacobi method"<<endl;
    cout << "Residual vector norm: " << endl;
    cout <<  residual(matrix, answer, vector, size) << endl;
    cout << endl;
    delete[]answer;
    for (int i = 0; i < size; i++) {
        delete[]C[i];
    }
    delete[]C;
    delete[]difference;
    delete[]temp;
    return 0;
}

int relaxationAndSeidelMethod(T parameter, const int size, T** matrix, T* vector) {
    T* answer = new T[size];
    T* temp = new T[size];
    for (int i = 0; i < size; i++) {
        answer[i] = 0;
        temp[i] = 0;
    }
    T sum1;
    T sum2;
    T A1=0; 
    T Ainf=0;
    T iter=0;
    T* difference = new T[size];
    for (int i = 0; i < size; i++) {
        difference[i] = 0;
    }
    do {
        iter++;
        for (int i = 0; i < size; i++) {
            sum1 = 0;
            sum2 = 0;
            for (int j = i + 1; j < size; j++) {
                sum2 += matrix[i][j] * temp[j];
            }

            for (int l = 0; l < i; l++) {
                sum1 += matrix[i][l] * answer[l];
            }
            if (parameter==1) {
                answer[i] = (vector[i] - sum1 - sum2) / matrix[i][i];
            }
            else {
                answer[i] = (1 - parameter) * temp[i] + parameter * (vector[i] - sum1 - sum2) / matrix[i][i];
            }
        }
        for (int k = 0; k < size; k++) {
            difference[k] = answer[k] - temp[k];
        }
        for (int k = 0; k < size; k++) {
            temp[k] = answer[k];
        }
    } while ((residual(matrix,answer,vector,size)/matrixNorm(size,matrix,A1,Ainf))>=epsilon);
    if (parameter == 1) {
        cout << "Solved by Seidel method: " << endl;
        for (int i = 0; i < size; i++) {
            cout << answer[i] << endl;
        }
        cout<<iter<<" iterations for Seidel method"<<endl;
    }
    else {
        cout << "Solved by relaxation method: " << endl;
        for (int i = 0; i < size; i++) {
            cout << answer[i] << endl;
        }
        cout<<iter<<" iterations for relaxation method"<<endl;
    }
    cout << "Residual vector norm: " << endl;
    cout <<  residual(matrix, answer, vector, size) << endl;
    cout << endl;
    delete[]answer;
    delete[]difference;
    delete[]temp;
    return 0;
}

int simpleIterationMethod(const int size, T** matrix, T* vector) {
    T tau = 0.0001;
    T iter = 0;
    T A1 = 0; T Ainf = 0;
    T normC = 0;
    T* answer = new T[size];
    T** C = new T * [size];
    for (int i = 0; i < size; i++) {
        C[i] = new T[size];
    }

    if (checkSimpleIterConvergence(size, matrix, C, tau, normC)) {
        cout << "Simple iterations method converges." << endl;
    }
    else {
        cout << "Matrix C norm is bigger than 1 for simple iterations method." << endl;
    }
    
    T* temp = new T[size];
    for (int i = 0; i < size; i++) {
        answer[i] = 0;
        temp[i] = 0;
    }
    T sum1;
    T* difference = new T[size];
    for (int i = 0; i < size; i++) {
        difference[i] = 0;
    }
    do {
        iter++;
        for (int i = 0; i < size; i++) {
            sum1 = 0;
            for (int j = 0; j < size; j++) {
                sum1 += matrix[i][j]*temp[j];
            }

            answer[i] = tau*(vector[i] - sum1) + temp[i];

        }
        for (int k = 0; k < size; k++) {
            difference[k] = answer[k] - temp[k];
        }
        for (int k = 0; k < size; k++) {
            temp[k] = answer[k];
        }
    } while ((residual(matrix,answer,vector,size)/matrixNorm(size,matrix,A1,Ainf))>=epsilon);
        cout << "Solved by simple iterations method: " << endl;
        for (int i = 0; i < size; i++) {
            cout << answer[i] << endl;
        }
    cout<<iter<<" iterations for simple iterations method"<<endl;
    cout << "Residual vector norm: " << endl;
    cout <<  residual(matrix, answer, vector, size) << endl;
    cout << endl;
    for (int i = 0; i < size; i++) {
        delete[]C[i];
    }
    delete[]C;
    delete[]answer;
    delete[]difference;
    delete[]temp;
    return 0;
}

//for 3-diagonal matrix
int relax3Method(const int size, T parameter, T* a,T* b, T* c, T* dvect) { 
    T* answer = new T[size];
    T* temp = new T[size];
    for (int i = 0; i < size; i++) {
        answer[i] = 0;
        temp[i] = 0;
    }
    T sum1;
    T sum2;
    T* difference = new T[size];
    for (int i = 0; i < size; i++) {
        difference[i] = 0;
    }
    do {
        answer[0]=parameter*(-c[0]*temp[1]+dvect[0])/b[0]+(1-parameter)*temp[0];
        for (int i = 1; i < size-1; i++) {
                sum2 = c[i] * temp[i+1];
                sum1 = a[i] * answer[i-1];
                answer[i] = (1 - parameter) * temp[i] + parameter * (dvect[i] - sum1 - sum2) / b[i];
        }
        answer[size-1]=(1-parameter)*temp[size-1]+parameter*(dvect[size-1]-a[size-1]*answer[size-2])/b[size-1];
        for (int k = 0; k < size; k++) {
            difference[k] = answer[k] - temp[k];
        }
        for (int k = 0; k < size; k++) {
            temp[k] = answer[k];
        }
    } while (vectorNorm(difference, size) >= epsilon);
        cout << "Solved by relaxation method: " << endl;
        for (int i = 0; i < size; i++) {
            cout << answer[i] << endl;
        }
    delete[]answer;
    delete[]difference;
    delete[]temp;
    return 0;
}


int main() {
    char* matrixFilename = new char[maxLength];
    char* vectorFilename = new char[maxLength];
    size_t sizet;
    readSettings("fileset.txt", sizet, matrixFilename, vectorFilename);
    const int size = sizet;
    T** matrix = new T* [size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new T[size];
    }
    T* vector = new T[size];
    T parameter = 0.1;
    
    
    readFileMatrixAndVector(matrixFilename, vectorFilename, matrix, vector, size);
    simpleIterationMethod(size, matrix, vector);
    if (checkConvergenceJacobiSeidel(size, matrix)) {
        cout << "Jacobi and relaxation methods converge." << endl;
        JacobiMethod(size, matrix, vector);
        relaxationAndSeidelMethod(1,size,matrix,vector);
        relaxationAndSeidelMethod(parameter, size, matrix, vector);
    }
    else {
        cout << "Jacobi and relaxation methods do not converge." << endl;
    }
    // const int sized = 200;
    // T* pdiag= new T[sized];
    // T* diag= new T[sized];
    // T* ndiag= new T[sized];
    // T* dvect= new T[sized];

//    for (int i=0;i<sized;i++){
//        diag[i]=8;
//        dvect[i]=i;
//        pdiag[i]=6;
//        ndiag[i]=1;
//    }
//    relax3Method(sized,parameter,pdiag,diag,ndiag,dvect);
    
//    for (int i=0;i<sized;i++){
//        diag[i]=8;
//        dvect[i]=i;
//        pdiag[i]=1;
//        ndiag[i]=6;
//    }
//    relax3Method(sized,parameter,pdiag,diag,ndiag,dvect);
    
    for (int i = 0; i < size; i++) {
        delete[]matrix[i];
    }
    delete[]matrix;
    delete[]vector;
    delete[]matrixFilename;
    delete[]vectorFilename;
    // delete[]diag;
    // delete[]ndiag;
    // delete[]pdiag;
    // delete[]dvect;
    return 0;
}

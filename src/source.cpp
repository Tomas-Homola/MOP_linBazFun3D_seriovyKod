#define _USE_MATH_DEFINES

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define R 6378.0
#define D 300.0
#define GM 398600.5
#define M 902
#define EPSILON 0.000000000001

using std::cos;
using std::sin;
using std::sqrt;

int main(int argc, char** argv)
{
    printf("Serial code\n\n");
    double* B = new double[M] {0.0};
    double* L = new double[M] {0.0};
    double Brad = 0.0, Lrad = 0.0, H = 0.0, u2n2 = 0.0;
    double temp = 0.0;

    //
    double gaussTable[7][4] = { 1.0 / 3.0 , 1.0 / 3.0 , 1.0 / 3.0 , 0.225,      // k = 1
                                0.79742699, 0.10128651, 0.10128651, 0.12593918, // k = 2
                                0.10128651, 0.79742699, 0.10128651, 0.12593918, // k = 3
                                0.10128651, 0.10128651, 0.79742699, 0.12593918, // k = 4
                                0.05971587, 0.47014206, 0.47014206, 0.13239415, // k = 5
                                0.47014206, 0.05971587, 0.47014206, 0.13239415, // k = 6
                                0.47014206, 0.47014206, 0.05971587, 0.13239415  // k = 7
    };

    // suradnice bodov X_i
    double* X = new double[M] {0.0};
    double* Y = new double[M] {0.0};
    double* Z = new double[M] {0.0};    

    // suradnicce normal v x_i
    double** n_x = new double*[M];
    double** n_y = new double*[M];
    double** n_z = new double*[M];

    for (int i = 0; i < M; i++)
    {
        n_x[i] = new double[6] {0.0};
        n_y[i] = new double[6] {0.0};
        n_z[i] = new double[6] {0.0};
    }

    // g vektor
    double* g = new double[M] {0.0};

    // Load GEOMETRY data
    printf("Loading geometry... ");
    FILE* file = nullptr;
    file = fopen("E:/_school/5_ZS/MOP/cv04_linearneBazoveFunkcie3D/BL-902.dat", "r");
    if (file == nullptr)
    {
        printf("Geometry file did not open\n");
        return -1;
    }
    for (int i = 0; i < M; i++)
    {
        int result = fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H, &g[i], &u2n2);
        //g[i] = -g[i] * 0.00001;
        g[i] = -GM / R;

        //g[i] = u2n2;
        Brad = B[i] * M_PI / 180.0;
        Lrad = L[i] * M_PI / 180.0;
        H = 0.0;
        X[i] = (R + H) * cos(Brad) * cos(Lrad);
        Y[i] = (R + H) * cos(Brad) * sin(Lrad);
        Z[i] = (R + H) * sin(Brad);

        //if (i < 10)
          //  printf("g[%d]: %.5lf\n", i, g[i]);
        //    printf("X[%d] = (%.2lf, %.2lf, %.2lf)\n", i, X_x[i], X_y[i], X_z[i]);
    }
    printf("done\n");
    fclose(file);   

    // Allocate space for array E
    int** E = new int*[M];
    for (int i = 0; i < M; i++)
    {
        // i - 1
        E[i] = new int[7] {0};
    }

    // Load ELEMENTS data
    printf("Loading elements data... ");
    file = fopen("E:/_school/5_ZS/MOP/cv04_linearneBazoveFunkcie3D/elem_902.dat", "r");
    for (int i = 0; i < M; i++)
    {
        int result = fscanf(file, "%d %d %d %d %d %d %d", &E[i][0], &E[i][1], &E[i][2], &E[i][3], &E[i][4], &E[i][5], &E[i][6]);
    }
    printf("done\n");
    fclose(file);

    /*for (int i = 0; i < 5; i++)
    {
        printf("%d: %d %d %d %d %d %d %d\n", i, E[i][0], E[i][1], E[i][2], E[i][3], E[i][4], E[i][5], E[i][6]);
    }*/

    // Compute areas of all triangles and normals to triangles
    int next2 = -1; // pomocna premenna pre vypocet indexu dalsieho suseda
    int next1 = -1;
    // vektory u a v -> u = E[j][t] - j, v = E[j][s] - j, w = u x v -> vektorovy sucin
    double u_x = 0.0;
    double u_y = 0.0;
    double u_z = 0.0;
    double v_x = 0.0;
    double v_y = 0.0;
    double v_z = 0.0;
    double w_x = 0.0;
    double w_y = 0.0;
    double w_z = 0.0;
    double w_norm = 0.0;
    double A_sum = 0.0;

    printf("Calculating triangles...");
    double** A = new double* [M];
    for (int j = 0; j < M; j++)
    {
        // kazdy j-ty support ma 4/6 trojuholnikov
        A[j] = new double[6] {0.0};
    
        // iteracia cez vsetky trojuholniky j-teho supportu -> ulozene v poli E na pozicii [j][0]
        for (int t = 1; t <= E[j][0]; t++) // od [1,6], riadky E su [0,7]
        {
            next1 = E[j][t] - 1;
            u_x = X[next1] - X[j];
            u_y = Y[next1] - Y[j];
            u_z = Z[next1] - Z[j];

            next2 = (t == E[j][0]) ? 1 : t + 1;

            v_x = X[E[j][next2] - 1] - X[j];
            v_y = Y[E[j][next2] - 1] - Y[j];
            v_z = Z[E[j][next2] - 1] - Z[j];

            w_x = u_y * v_z - v_y * u_z;
            w_y = u_z * v_x - v_z * u_x;
            w_z = u_x * v_y - v_x * u_y;

            // velkost vektora w -> ...
            w_norm = std::sqrt(w_x * w_x + w_y * w_y + w_z * w_z);

            // ... polovica je plocha trojuholnika A[j][t-1]
            temp = w_norm / 2.0;
            A[j][t - 1] = temp;
            A_sum += temp;

            // normala ku trojuholniku je normovany vektor w
            n_x[j][t - 1] = w_x / w_norm;
            n_y[j][t - 1] = w_y / w_norm;
            n_z[j][t - 1] = w_z / w_norm;
        }
    }
    printf("done\n\n");
    double A_earth = 4.0 * M_PI * R * R;
    printf("Sphere area: %lf km^2\nrelative error: %lf\n", A_sum / 3.0, abs(A_sum / 3.0 - A_earth) / A_earth);
    printf("Sphere/Earth: %lf\n\n", (A_sum / 3.0) / A_earth);

    printf("Computing gauss points... ");
    // suradnice gaussovych bodov
    double*** x = new double** [M];
    double*** y = new double** [M];
    double*** z = new double** [M];

    for (int j = 0; j < M; j++)
    {
        x[j] = new double* [6];
        y[j] = new double* [6];
        z[j] = new double* [6];

        for (int t = 1; t <= 6; t++)
        {
            x[j][t - 1] = new double[7] {0.0};
            y[j][t - 1] = new double[7] {0.0};
            z[j][t - 1] = new double[7] {0.0};

            for (int k = 0; k < 7; k++)
            {
                next1 = E[j][t] - 1;
                next2 = (t == E[j][0]) ? 1 : t + 1;

                x[j][t - 1][k] = X[j] * gaussTable[k][0] + X[next1] * gaussTable[k][1] + X[next2] * gaussTable[k][2];
                y[j][t - 1][k] = Y[j] * gaussTable[k][0] + Y[next1] * gaussTable[k][1] + Y[next2] * gaussTable[k][2];
                z[j][t - 1][k] = Z[j] * gaussTable[k][0] + Z[next1] * gaussTable[k][1] + Z[next2] * gaussTable[k][2];
            }
        }
    }

    printf("done\n");
    //for (int i = 0; i < N; i++) // set constant g values
    //    g[i] = -(GM) / (R * R);

    // vytvorenie matice systemu rovnic
    /*double* A = new double[N * N] {0.0};*/
    int ij = -1;





    //########## BCGS linear solver ##########//

    //double* sol = new double[N]; // vektor x^0 -> na ukladanie riesenia systemu
    //double* r_hat = new double[N]; // vektor \tilde{r} = b - A.x^0;
    //double* r = new double[N]; // vektor pre rezidua
    //double* p = new double[N]; // pomocny vektor na update riesenia
    //double* v = new double[N]; // pomocny vektor na update riesenia
    //double* s = new double[N]; // pomocny vektor na update riesenia
    //double* t = new double[N]; // pomocny vektor na update riesenia

    //double beta = 0.0;
    //double rhoNew = 1.0;
    //double rhoOld = 0.0;
    //double alpha = 1.0;
    //double omega = 1.0;

    //double tempDot = 0.0;
    //double tempDot2 = 0.0;
    //double sNorm = 0.0;

    //int MAX_ITER = 1000;
    //double TOL = 1.0E-6;
    //int iter = 1;

    //double rezNorm = 0.0;
    //for (int i = 0; i < N; i++) // set all to zero
    //{
    //    sol[i] = 0.0;
    //    p[i] = 0.0; // = 0
    //    v[i] = 0.0; // = 0
    //    s[i] = 0.0;
    //    t[i] = 0.0;

    //    r[i] = g[i];
    //    r_hat[i] = g[i];
    //    rezNorm += r[i] * r[i];

    //}

    //printf("||r0||: %.10lf\n", sqrt(rezNorm));
    //rezNorm = 0.0;

    //do
    //{
    //    rhoOld = rhoNew; // save previous rho_{i-2}
    //    rhoNew = 0.0; // compute new rho_{i-1}
    //    for (int i = 0; i < N; i++) // dot(r_hat, r)
    //        rhoNew += r_hat[i] * r[i];

    //    if (rhoNew == 0.0)
    //        return -1;

    //    if (iter == 1)
    //    {
    //        //printf("iter 1 setup\n");
    //        for (int i = 0; i < N; i++)
    //            p[i] = r[i];
    //    }
    //    else
    //    {
    //    beta = (rhoNew / rhoOld) * (alpha / omega);
    //    for (int i = 0; i < N; i++) // update vector p^(i)
    //        p[i] = r[i] + beta * (p[i] - omega * v[i]);
    //    }

    //    // compute vector v = A.p
    //    for (int i = 0; i < N; i++)
    //    {
    //        v[i] = 0.0;
    //        for (int j = 0; j < N; j++)
    //        {
    //            ij = i * N + j;
    //            v[i] += A[ij] * p[j];
    //        }
    //    }

    //    // compute alpha
    //    tempDot = 0.0;
    //    for (int i = 0; i < N; i++)
    //        tempDot += r_hat[i] * v[i];

    //    alpha = rhoNew / tempDot;

    //    // compute vektor s
    //    sNorm = 0.0;
    //    for (int i = 0; i < N; i++)
    //    {
    //        s[i] = r[i] - alpha * v[i];
    //        sNorm += s[i] * s[i];
    //    }

    //    sNorm = sqrt(sNorm);
    //    if (sNorm < TOL) // check if ||s|| is small enough
    //    {
    //        for (int i = 0; i < N; i++) // update solution x
    //            sol[i] = sol[i] + alpha * p[i];

    //        printf("BCGS stop:   ||s||(= %.10lf) is small enough, iter: %3d\n", sNorm, iter);
    //        break;
    //    }

    //    // compute vector t = A.s
    //    for (int i = 0; i < N; i++)
    //    {
    //        t[i] = 0.0;
    //        for (int j = 0; j < N; j++)
    //        {
    //            ij = i * N + j;
    //            t[i] += A[ij] * s[j];
    //        }
    //    }

    //    // compute omega
    //    tempDot = 0.0; tempDot2 = 0.0;
    //    for (int i = 0; i < N; i++)
    //    {
    //        tempDot += t[i] * s[i];
    //        tempDot2 += t[i] * t[i];
    //    }
    //    omega = tempDot / tempDot2;

    //    rezNorm = 0.0;
    //    for (int i = 0; i < N; i++)
    //    {
    //        sol[i] = sol[i] + alpha * p[i] + omega * s[i]; // update solution x
    //        r[i] = s[i] - omega * t[i]; // compute new residuum vector
    //        rezNorm += r[i] * r[i]; // compute residuum norm
    //    }

    //    rezNorm = sqrt(rezNorm);
    //    printf("iter: %3d    ||r||: %.10lf\n", iter, rezNorm);

    //    if (rezNorm < TOL)
    //    {
    //        printf("BCGS stop iter: ||r|| is small enough\n");
    //        break;
    //    }

    //    iter++;

    //} while ((iter < MAX_ITER) && (rezNorm > TOL));

    //delete[] r_hat;
    //delete[] r;
    //delete[] p;
    //delete[] v;
    //delete[] s;
    //delete[] t;

    ////########## EXPORT DATA ##########//
    //file = fopen("../outCorrect_serial.dat", "w");
    //if (file == nullptr)
    //{
    //    printf("data export failed\n");
    //    return -1;
    //}

    //printf("solution export started... ");
    //double ui = 0.0, Gij = 0.0;
    //for (int i = 0; i < N; i++)
    //{
    //    ui = 0.0;
    //    for (int j = 0; j < N; j++) // compute solution u(X_i)
    //    {
    //        r_x = X_x[i] - s_x[j];
    //        r_y = X_y[i] - s_y[j];
    //        r_z = X_z[i] - s_z[j];

    //        rNorm = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

    //        Gij = 1.0 / (4.0 * M_PI * rNorm);

    //        ui += sol[j] * Gij;
    //    }

    //    fprintf(file, "%.5lf\t%.5lf\t%.5lf\n", B[i], L[i], ui);
    //}

    //fclose(file);
    //printf("done\n");

    delete[] X;
    delete[] Y;
    delete[] Z;
    delete[] n_x;
    delete[] n_y;
    delete[] n_z;
    delete[] g;
}
#define _USE_MATH_DEFINES

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <memory>
#include <chrono>

#define R 6378.0
#define GM 398600.5
//#define M 902
#define M 3602
//#define M 160002
#define EPSILON 0.000000000001
#define TOL 1.0E-6
#define MAX_ITER  1000

#define pos_wk 3
#define pos_eta1k 0

using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::log;
using std::sqrt;

double Gii_ln(double beta_t, double theta_t)
{
    double result = 0.0;
    double nominator = tan((beta_t + theta_t) / 2.0);
    double denominator = tan(beta_t / 2.0);

    result = log(nominator / denominator);
    return result;
}

double angle(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double angle = 0.0;

    double dotProd = x1 * x2 + y1 * y2 + z1 * z2;
    double norm1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double norm2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

    angle = acos(dotProd / (norm1 * norm2));

    return angle;
}

int main(int argc, char** argv)
{
    int nprocs = 6;
    omp_set_num_threads(nprocs);

    auto start = std::chrono::high_resolution_clock::now();

    if (nprocs == 1)
        printf("Serial code\n\n");
    else
        printf("Parallel code (nprocs: %d)\n\n", nprocs);

    double* B = new double[M] {0.0};
    double* L = new double[M] {0.0};
    double Brad = 0.0, Lrad = 0.0, H = 0.0, u2n2 = 0.0;
    double temp = 0.0;
    int i = 0; int j = 0; int t = 0; int k = 0;

    // hodnoty pre gaussove body [eta1k, eta2k, eta3k, wk]
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
    double** n_x = new double* [M];
    double** n_y = new double* [M];
    double** n_z = new double* [M];

    for (i = 0; i < M; i++)
    {
        n_x[i] = new double[6] {0.0};
        n_y[i] = new double[6] {0.0};
        n_z[i] = new double[6] {0.0};
    }

    // q vektor
    double* q = new double[M] {0.0};

    // Load GEOMETRY data
    printf("Loading geometry... ");
    FILE* file = nullptr;
    //file = fopen("BL-902.dat", "r");
    file = fopen("BL-3602.dat", "r");
    //file = fopen("BL-160002.dat", "r");

    if (file == nullptr)
    {
        printf("Geometry file did not open\n");
        return -1;
    }
    for (i = 0; i < M; i++)
    {
        int result = fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H, &q[i], &u2n2);
        //q[i] = q[i] * 0.00001;
        q[i] = GM / (R * R);

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
    int** E = new int* [M];
    for (i = 0; i < M; i++)
    {
        // i - 1
        E[i] = new int[7] {0};
    }

    // Load ELEMENTS data
    printf("Loading elements data... ");
    //file = fopen("elem_902.dat", "r");
    file = fopen("elem_3602.dat", "r");
    //file = fopen("elem_160002.dat", "r");

    for (i = 0; i < M; i++)
    {
        int result = fscanf(file, "%d %d %d %d %d %d %d", &E[i][0], &E[i][1], &E[i][2], &E[i][3], &E[i][4], &E[i][5], &E[i][6]);
    }
    printf("done\n");
    fclose(file);

    // Compute areas of all triangles and normals to triangles
    int next2 = -1; // pomocna premenna pre vypocet indexu dalsieho suseda
    int s = -1;
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
    double Asum = 0.0;
    double tSum = 0.0;
    double beta_t = 0.0;
    double theta_t = 0.0;
    double l_t = 0.0;

    printf("Calculating triangles...");
    double** A = new double* [M];
    double* Gdiag = new double[M] {0.0};
#pragma omp parallel for private(t, next1, u_x, u_y, u_z, s, next2, v_x, v_y, v_z, w_x, w_y, w_z, w_norm, temp, Asum, theta_t, beta_t, l_t, tSum)
    for (j = 0; j < M; j++)
    {
        // kazdy j-ty support ma 4/6 trojuholnikov
        A[j] = new double[6] {0.0};

        tSum = 0.0;
        // iteracia cez vsetky trojuholniky j-teho supportu -> ulozene v poli E na pozicii [j][0]
        for (t = 1; t <= E[j][0]; t++) // od [1,6], riadky E su [0,7]
        {
            next1 = E[j][t] - 1;
            u_x = X[next1] - X[j];
            u_y = Y[next1] - Y[j];
            u_z = Z[next1] - Z[j];

            s = (t == E[j][0]) ? 1 : t + 1;
            next2 = E[j][s] - 1;

            v_x = X[next2] - X[j];
            v_y = Y[next2] - Y[j];
            v_z = Z[next2] - Z[j];

            w_x = u_y * v_z - v_y * u_z;
            w_y = u_z * v_x - v_z * u_x;
            w_z = u_x * v_y - v_x * u_y;

            // velkost vektora w -> ...
            w_norm = std::sqrt(w_x * w_x + w_y * w_y + w_z * w_z);

            // ... polovica je plocha trojuholnika A[j][t-1]
            temp = w_norm / 2.0;
            A[j][t - 1] = temp;
            //Asum += temp;

            // normala ku trojuholniku je normovany vektor w
            n_x[j][t - 1] = -w_x / w_norm;
            n_y[j][t - 1] = -w_y / w_norm;
            n_z[j][t - 1] = -w_z / w_norm;

            // VYPOCET Gii
            // uhol theta_t
            theta_t = angle(u_x, u_y, u_z, v_x, v_y, v_z);

            w_x = X[next2] - X[next1];
            w_y = Y[next2] - Y[next1];
            w_z = Z[next2] - Z[next1];

            // uhol beta_t
            beta_t = angle(-u_x, -u_y, -u_z, w_x, w_y, w_z);

            // dlzka l_t
            l_t = sqrt(w_x * w_x + w_y * w_y + w_z * w_z);

            temp = Gii_ln(beta_t, theta_t);
            tSum += (A[j][t - 1] / l_t) * temp;
        } // END TRIANGLES

        temp = tSum / (4.0 * M_PI);
        Gdiag[j] = temp;
    }

    printf("done\n\n");
    /*double A_earth = 4.0 * M_PI * R * R;
    printf("Sphere area: %.10lf km^2\n", A_sum / 3.0);
    printf("Sphere/Earth: %lf\n\n", (A_sum / 3.0) / A_earth);*/

    //for (int i = 0; i < 7; i++)
    //{
    //    printf("moje Gdiag[%d]: %.5lf\n", i, Gdiag[i]);
    //}

    printf("Computing gauss points... ");
    // suradnice gaussovych bodov
    double*** x = new double** [M];
    double*** y = new double** [M];
    double*** z = new double** [M];

#pragma omp parallel for private(t, k, next1, s, next2)
    for (j = 0; j < M; j++)
    {
        x[j] = new double* [6] {nullptr};
        y[j] = new double* [6] {nullptr};
        z[j] = new double* [6] {nullptr};

        for (t = 1; t <= E[j][0]; t++)
        {
            x[j][t - 1] = new double[7] {0.0};
            y[j][t - 1] = new double[7] {0.0};
            z[j][t - 1] = new double[7] {0.0};

            for (k = 0; k < 7; k++)
            {
                next1 = E[j][t] - 1;
                s = (t == E[j][0]) ? 1 : t + 1;
                next2 = E[j][s] - 1;

                x[j][t - 1][k] = X[j] * gaussTable[k][0] + X[next1] * gaussTable[k][1] + X[next2] * gaussTable[k][2];
                y[j][t - 1][k] = Y[j] * gaussTable[k][0] + Y[next1] * gaussTable[k][1] + Y[next2] * gaussTable[k][2];
                z[j][t - 1][k] = Z[j] * gaussTable[k][0] + Z[next1] * gaussTable[k][1] + Z[next2] * gaussTable[k][2];
            }
        }
    }

    printf("done\n\n");

    // compute matrix F
    double F_ij = 0.0;
    double G_ij = 0.0;
    double tSum_F = 0.0;
    double kSum_F = 0.0;
    double tSum_G = 0.0;
    double kSum_G = 0.0;
    double r_ijk = 0.0;
    double xDist = 0.0; double yDist = 0.0; double zDist = 0.0;
    double w_k = 0.0; double psi_k = 0.0;
    double K_tij = 0.0;
    double r_x = 0.0; double r_y = 0.0; double r_z = 0.0;

    printf("Computing matrix F... ");
    double** F = new double* [M];
    //double* G = new double[M] {0.0};
    double* rhs = new double[M] {0.0};

    for (i = 0; i < M; i++)
    {
        F[i] = new double[M] {0.0};

        F[i][i] = 1.0;
    }

#pragma omp parallel for private(j, t, k, F_ij, G_ij, tSum_F, tSum_G, xDist, yDist, zDist, r_ijk, w_k, psi_k, temp, kSum_F, kSum_G, K_tij, r_x, r_y, r_z)
    for (i = 0; i < M; i++)
    {
        F_ij = 0.0;
        for (j = 0; j < M; j++) // iteracie cez vsetky supporty "j"
        {
            if (i != j) // iba nesingularne
            {
                tSum_F = 0.0;
                tSum_G = 0.0;
                for (t = 1; t <= E[j][0]; t++) // iteracie cez vsetky trojuholniky supportu "j"
                {
                    kSum_F = 0.0;
                    kSum_G = 0.0;
                    for (k = 0; k < 7; k++) // iteracie cez vsetky Gaussove body
                    {
                        // x[j][t - 1][k]
                        xDist = x[j][t - 1][k] - X[i];
                        yDist = y[j][t - 1][k] - Y[i];
                        zDist = z[j][t - 1][k] - Z[i];

                        r_ijk = sqrt(xDist * xDist + yDist * yDist + zDist * zDist);
                        w_k = gaussTable[k][pos_wk];
                        psi_k = gaussTable[k][pos_eta1k];

                        temp = w_k * psi_k / (r_ijk * r_ijk * r_ijk);
                        kSum_F += temp;

                        temp = w_k * psi_k / r_ijk;
                        kSum_G += temp;
                    } // END GAUSS

                    // n_x[j][t - 1] = w_x / w_norm;
                    // vektor z i-teho do j-teho kolokacneho bodu
                    r_x = X[i] - X[j];
                    r_y = Y[i] - Y[j];
                    r_z = Z[i] - Z[j];

                    // skalarny sucin medzi vektorom r_ij a normalou n^t_j
                    K_tij = r_x * n_x[j][t - 1] + r_y * n_y[j][t - 1] + r_z * n_z[j][t - 1];

                    tSum_F += A[j][t - 1] * K_tij * kSum_F;

                    tSum_G += A[j][t - 1] * kSum_G;
                } // END TRIANGLES

                F_ij = tSum_F / (4.0 * M_PI);

                F[i][j] = F_ij;
                F[i][i] -= F_ij;

                G_ij = tSum_G / (4.0 * M_PI);

                rhs[i] += G_ij * q[j];

            } // END IF NON-SINGULAR
            else
            {
                rhs[i] += Gdiag[j] * q[j];
            }
        }
    }
    printf("done\n\n");

    /*for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            printf("%.4lf\t", F[i][j]);
        }
        printf("\n");
    }*/

    //for (int i = 0; i < 9; i++)
    //{
    //    //printf("G[%d][%d]: %.5lf\n", i, i, Gdiag[i]);
    //    printf("rhs[%d]: %.5lf\n",i, rhs[i]);

    //}

    //exit(-2);
    //########## BCGS linear solver ##########//

    double* sol = new double[M]; // vektor x^0 -> na ukladanie riesenia systemu
    double* r_hat = new double[M]; // vektor \tilde{r} = b - A.x^0;
    double* r = new double[M]; // vektor pre rezidua
    double* p = new double[M]; // pomocny vektor na update riesenia
    double* v = new double[M]; // pomocny vektor na update riesenia
    double* sv = new double[M]; // pomocny vektor na update riesenia
    double* tv = new double[M]; // pomocny vektor na update riesenia

    double beta = 0.0;
    double rhoNew = 1.0;
    double rhoOld = 0.0;
    double alpha = 1.0;
    double omega = 1.0;

    double tempDot = 0.0;
    double tempDot2 = 0.0;
    double sNorm = 0.0;

    int iter = 1;

    double rezNorm = 0.0;
    for (i = 0; i < M; i++) // set all to zero
    {
        sol[i] = 0.0;
        p[i] = 0.0; // = 0
        v[i] = 0.0; // = 0
        sv[i] = 0.0;
        tv[i] = 0.0;

        r[i] = rhs[i];
        r_hat[i] = rhs[i];
        rezNorm += r[i] * r[i];

    }

    printf("||r0||: %.10lf\n", sqrt(rezNorm));
    rezNorm = 0.0;

    do
    {
        rhoOld = rhoNew; // save previous rho_{i-2}
        rhoNew = 0.0; // compute new rho_{i-1}
        for (i = 0; i < M; i++) // dot(r_hat, r)
            rhoNew += r_hat[i] * r[i];

        if (rhoNew == 0.0)
            return -1;

        if (iter == 1)
        {
            //printf("iter 1 setup\n");
            for (i = 0; i < M; i++)
                p[i] = r[i];
        }
        else
        {
            beta = (rhoNew / rhoOld) * (alpha / omega);
            for (i = 0; i < M; i++) // update vector p^(i)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // compute vector v = A.p
#pragma omp parallel for private(j)
        for (i = 0; i < M; i++)
        {
            v[i] = 0.0;
            for (j = 0; j < M; j++)
            {
                v[i] += F[i][j] * p[j];
            }
        }

        // compute alpha
        tempDot = 0.0;
        for (i = 0; i < M; i++)
            tempDot += r_hat[i] * v[i];

        alpha = rhoNew / tempDot;

        // compute vektor s
        sNorm = 0.0;
        for (i = 0; i < M; i++)
        {
            sv[i] = r[i] - alpha * v[i];
            sNorm += sv[i] * sv[i];
        }

        sNorm = sqrt(sNorm);
        if (sNorm < TOL) // check if ||s|| is small enough
        {
            for (i = 0; i < M; i++) // update solution x
                sol[i] = sol[i] + alpha * p[i];

            printf("BCGS stop:   ||s||(= %.10lf) is small enough, iter: %3d\n", sNorm, iter);
            break;
        }

        // compute vector t = A.s
#pragma omp parallel for private(j)
        for (i = 0; i < M; i++)
        {
            tv[i] = 0.0;
            for (j = 0; j < M; j++)
            {
                tv[i] += F[i][j] * sv[j];
            }
        }

        // compute omega
        tempDot = 0.0; tempDot2 = 0.0;
        for (i = 0; i < M; i++)
        {
            tempDot += tv[i] * sv[i];
            tempDot2 += tv[i] * tv[i];
        }
        omega = tempDot / tempDot2;

        rezNorm = 0.0;
        for (i = 0; i < M; i++)
        {
            sol[i] = sol[i] + alpha * p[i] + omega * sv[i]; // update solution x
            r[i] = sv[i] - omega * tv[i]; // compute new residuum vector
            rezNorm += r[i] * r[i]; // compute residuum norm
        }

        rezNorm = sqrt(rezNorm);
        printf("iter: %3d    ||r||: %.10lf\n", iter, rezNorm);

        if (rezNorm < TOL)
        {
            printf("BCGS stop iter: ||r|| is small enough\n");
            break;
        }

        iter++;

    } while ((iter < MAX_ITER) && (rezNorm > TOL));

    delete[] r_hat;
    delete[] r;
    delete[] p;
    delete[] v;
    delete[] sv;
    delete[] tv;

    double residuum = 0.0;
    //printf("\n\nGMR: %.4lf\n\n", GM / R);
    for (i = 0; i < M; i++)
    {
        residuum += (sol[i] - (GM / R)) * (sol[i] - (GM / R));
    
        if (i < 6)
            printf("sol[%d]: %.4lf\tGM/R: %.4lf\n", i, sol[i], GM / R);
    }

    printf("\n902 res: %.8lf\n", sqrt(residuum / M));
    exit(1);

    //########## EXPORT DATA ##########//
    file = fopen("sol_.dat", "w");
    if (file == nullptr)
    {
        printf("data export failed\n");
        return -1;
    }

    printf("solution export started... ");
    for (i = 0; i < M; i++)
    {
        fprintf(file, "%.6lf\t%.6lf\t%.6lf\n", B[i], L[i], sol[i]);
    }

    fclose(file);
    printf("done\n");

    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    printf("Duration %.2lf s\n", (double)time.count() / 1000.0);

    // Memory clean-up
    delete[] sol;
    delete[] B;
    delete[] L;
    delete[] X;
    delete[] Y;
    delete[] Z;
    for (i = 0; i < M; i++)
    {
        delete[] n_x[i];
        delete[] n_y[i];
        delete[] n_z[i];
    }
    delete[] n_x;
    delete[] n_y;
    delete[] n_z;
    delete[] q;
    for (i = 0; i < M; i++)
        delete[] E[i];
    delete[] E;
    delete[] Gdiag;
    for (j = 0; j < M; j++)
        delete[] A[j];
    delete[] A;

    for (j = 0; j < M; j++)
    {
        for (t = 0; t < 6; t++)
        {
            delete[] x[j][t];
            delete[] y[j][t];
            delete[] z[j][t];
        }
        delete[] x[j];
        delete[] y[j];
        delete[] z[j];
    }
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] rhs;
    for (i = 0; i < M; i++)
        delete[] F[i];
    delete[] F;
}
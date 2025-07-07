using MathNet.Numerics.LinearAlgebra;

namespace PracticeLibrary
{
    public class Solve
    {
        public static double[,,] SolveGeneral(int nx, int ny, int nz, double lx, double ly, double lz, double tMax, double dt, double alpha, Func<double, double, double, double> initialCondition, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients)
        {
            double dx = lx / (nx - 1);
            double dy = ly / (ny - 1);
            double dz = lz / (nz - 1);

            double rx = alpha * dt / (2.0 * dx * dx);
            double ry = alpha * dt / (2.0 * dy * dy);
            double rz = alpha * dt / (2.0 * dz * dz);

            double[,,] solution = new double[nx, ny, nz];
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nz; k++)
                    {
                        double x = i * dx;
                        double y = j * dy;
                        double z = k * dz;
                        solution[i, j, k] = initialCondition(x, y, z);
                    }
                }
            }

            Console.WriteLine($"Проверка начальных условий:");
            for (int k = 0; k < nz; k++)
                Console.WriteLine($"Слой z={k}: {solution[0, 0, k]} {solution[1, 0, k]} {solution[2, 0, k]}");

            double totalMass = 0;
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                        totalMass += solution[i, j, k] * dx * dy * dz;
            Console.WriteLine($"Total mass: {totalMass}\tt = 0");

            for (double t = 0; t < tMax; t += dt)
            {
                SolveX(solution, alpha, dt, rx, nx, ny, nz, dx, dy, dz, boundaryConditions, robinCoefficients, t);
                SolveY(solution, alpha, dt, ry, nx, ny, nz, dx, dy, dz, boundaryConditions, robinCoefficients, t);
                SolveZ(solution, alpha, dt, rz, nx, ny, nz, dx, dy, dz, boundaryConditions, robinCoefficients, t);

                Console.WriteLine($"После шага t={t}:");
                for (int k = 0; k < nz; k++)
                    Console.WriteLine($"Слой z={k}: {solution[0, 0, k]:F6} {solution[1, 0, k]:F6} {solution[2, 0, k]:F6}");

                totalMass = 0;
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k < nz; k++)
                            totalMass += solution[i, j, k] * dx * dy * dz;
                Console.WriteLine($"Total mass: {totalMass}\tt = {t}");
            }

            return solution;
        }

        private static void SolveX(double[,,] solution, double alpha, double dt, double rx, int nx, int ny, int nz, double dx, double dy, double dz, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients, double t)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nz; k++)
                {
                    Matrix<double> matrix = Create.TridiagonalMatrixX(nx, rx);
                    Vector<double> vector = Create.RightHandSideX(solution, alpha, dt, rx, nx, j, k, dx, dy, dz, t);


                    Apply.BoundaryConditionsX(matrix, vector, rx, nx, j, dx, boundaryConditions, robinCoefficients);

                    Vector<double> x = matrix.Solve(vector);

                    for (int i = 0; i < nx; i++)
                    {
                        solution[i, j, k] = x[i];
                    }
                }
            }
        }

        private static void SolveY(double[,,] solution, double alpha, double dt, double ry, int nx, int ny, int nz, double dx, double dy, double dz, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients, double t)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int k = 0; k < nz; k++)
                {
                    Matrix<double> matrix = Create.TridiagonalMatrixY(ny, ry);
                    Vector<double> vector = Create.RightHandSideY(solution, alpha, dt, ry, ny, i, k, dx, dy, dz, t);


                    Apply.BoundaryConditionsY(matrix, vector, ry, ny, dy, boundaryConditions, robinCoefficients);

                    Vector<double> y = matrix.Solve(vector);

                    for (int j = 0; j < ny; j++)
                    {
                        solution[i, j, k] = y[j];
                    }
                }
            }
        }

        private static void SolveZ(double[,,] solution, double alpha, double dt, double rz, int nx, int ny, int nz, double dx, double dy, double dz, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients, double t)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    Matrix<double> matrix = Create.TridiagonalMatrixZ(nz, rz);
                    Vector<double> vector = Create.RightHandSideZ(solution, alpha, dt, rz, nz, i, j, dx, dy, dz, t);


                    Apply.BoundaryConditionsZ(matrix, vector, rz, nz, dz, boundaryConditions, robinCoefficients);

                    Vector<double> z = matrix.Solve(vector);

                    for (int k = 0; k < nz; k++)
                    {
                        solution[i, j, k] = z[k];
                    }
                }
            }
        }
    }
}

using MathNet.Numerics.LinearAlgebra;

namespace PracticeLibrary
{
    public class Solve_2
    {
        private static double CalculateTotalMass(double[,,] solution, double dx, double dy, double dz)
        {
            double totalMass = 0;
            for (int i = 0; i < solution.GetLength(0); i++)
                for (int j = 0; j < solution.GetLength(1); j++)
                    for (int k = 0; k < solution.GetLength(2); k++)
                        totalMass += solution[i, j, k] * dx * dy * dz;
            return totalMass;
        }
        public static (double[,,] algaeSolution, double[,,] shrimpsSolution) SolveGeneral(
            int nx,
            int ny,
            int nz,
            double lx,
            double ly,
            double lz,
            double tMax,
            double dt,
            double algaeDiffusion,
            double shrimpDiffusion,
            Func<double, double, double, double> initialAlgae,
            Func<double, double, double, double> initialShrimps,
            Dictionary<string, (int type, double value)> algaeBoundaryConditions,
            Dictionary<string, (int type, double value)> shrimpsBoundaryConditions,
            Dictionary<string, (double alpha, double beta)> robinCoefficients)
        {
            double algaeGrowthRate = 0;    // Скорость роста водорослей
            double shrimpDeathRate = 0;    // Смертность креветок
            double interactionRate = 1;    // Скорость поедания водорослей креветками
            double conversionRate = 1;    // Эффективность преобразования водорослей в креветок

            double dx = lx / (nx - 1);
            double dy = ly / (ny - 1);
            double dz = lz / (nz - 1);

            double rx_aglae = algaeDiffusion * dt / (2.0 * dx * dx);
            double ry_aglae = algaeDiffusion * dt / (2.0 * dy * dy);
            double rz_aglae = algaeDiffusion * dt / (2.0 * dz * dz);

            // Для креветок
            double rx_shrimp = shrimpDiffusion * dt / (2.0 * dx * dx);
            double ry_shrimp = shrimpDiffusion * dt / (2.0 * dy * dy);
            double rz_shrimp = shrimpDiffusion * dt / (2.0 * dz * dz);

            double[,,] algaeSolution = new double[nx, ny, nz];
            double[,,] shrimpsSolution = new double[nx, ny, nz];
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nz; k++)
                    {
                        double x = i * dx;
                        double y = j * dy;
                        double z = k * dz;
                        algaeSolution[i, j, k] = initialAlgae(x, y, z);
                        shrimpsSolution[i, j, k] = initialShrimps(x, y, z);

                    }
                }
            }

            // Вывод информации о начальных условиях
            Console.WriteLine("Начальные условия:");
            Console.WriteLine($"Общая масса водорослей: {CalculateTotalMass(algaeSolution, dx, dy, dz)}");
            Console.WriteLine($"Общая масса креветок: {CalculateTotalMass(shrimpsSolution, dx, dy, dz)}");
            Console.WriteLine($"Общая масса: {CalculateTotalMass(shrimpsSolution, dx, dy, dz) + CalculateTotalMass(algaeSolution, dx, dy, dz)}");


            for (double t = 0; t < tMax; t += dt)
            {
                SolveX(algaeSolution, algaeDiffusion, dt, rx_aglae, nx, ny, nz, dx, dy, dz,
                    algaeBoundaryConditions, robinCoefficients, t, shrimpsSolution,
                    (a, s) => a * algaeGrowthRate * (1 - a) - interactionRate * a * s);
                SolveY(algaeSolution, algaeDiffusion, dt, ry_aglae, nx, ny, nz, dx, dy, dz,
                    algaeBoundaryConditions, robinCoefficients, t, shrimpsSolution,
                    (a, s) => a * algaeGrowthRate * (1 - a) - interactionRate * a * s);
                SolveZ(algaeSolution, algaeDiffusion, dt, rz_aglae, nx, ny, nz, dx, dy, dz,
                    algaeBoundaryConditions, robinCoefficients, t, shrimpsSolution,
                    (a, s) => a * algaeGrowthRate * (1 - a) - interactionRate * a * s);

                SolveX(shrimpsSolution, shrimpDiffusion, dt, rx_shrimp, nx, ny, nz, dx, dy, dz,
              shrimpsBoundaryConditions, robinCoefficients, t, algaeSolution,
              (s, a) => -shrimpDeathRate * s + conversionRate * a * s);

                SolveY(shrimpsSolution, shrimpDiffusion, dt, ry_shrimp, nx, ny, nz, dx, dy, dz,
                      shrimpsBoundaryConditions, robinCoefficients, t, algaeSolution,
                      (s, a) => -shrimpDeathRate * s + conversionRate * a * s);

                SolveZ(shrimpsSolution, shrimpDiffusion, dt, rz_shrimp, nx, ny, nz, dx, dy, dz,
                      shrimpsBoundaryConditions, robinCoefficients, t, algaeSolution,
                      (s, a) => -shrimpDeathRate * s + conversionRate * a * s);

                Console.WriteLine($"После шага t={t}:");
                Console.WriteLine($"Общая масса водорослей: {CalculateTotalMass(algaeSolution, dx, dy, dz)}");
                Console.WriteLine($"Общая масса креветок: {CalculateTotalMass(shrimpsSolution, dx, dy, dz)}");
                Console.WriteLine($"Общая масса: {CalculateTotalMass(shrimpsSolution, dx, dy, dz) + CalculateTotalMass(algaeSolution, dx, dy, dz)}");
            }

            return (algaeSolution, shrimpsSolution);
        }

        private static void SolveX(
        double[,,] solution,
        double alpha, double dt, double rx,
        int nx, int ny, int nz,
        double dx, double dy, double dz,
        Dictionary<string, (int type, double value)> boundaryConditions,
        Dictionary<string, (double alpha, double beta)> robinCoefficients,
        double t,
        double[,,] otherSolution,
        Func<double, double, double> interactionTerm)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nz; k++)
                {
                    Matrix<double> matrix = Create.TridiagonalMatrixX(nx, rx);
                    Vector<double> vector = Create.RightHandSideX(solution, alpha, dt, rx, nx, j, k, dx, dy, dz, t);

                    // Добавляем член взаимодействия
                    for (int i = 0; i < nx; i++)
                    {
                        vector[i] += dt * interactionTerm(solution[i, j, k], otherSolution[i, j, k]);
                    }

                    Apply.BoundaryConditionsX(matrix, vector, rx, nx, j, dx, boundaryConditions, robinCoefficients);
                    Vector<double> x = matrix.Solve(vector);

                    for (int i = 0; i < nx; i++)
                    {
                        solution[i, j, k] = Math.Max(0, x[i]); // Гарантируем неотрицательность
                    }
                }
            }

        }

        private static void SolveY(
        double[,,] solution,
        double alpha, double dt, double ry,
        int nx, int ny, int nz,
        double dx, double dy, double dz,
        Dictionary<string, (int type, double value)> boundaryConditions,
        Dictionary<string, (double alpha, double beta)> robinCoefficients,
        double t,
        double[,,] otherSolution,
        Func<double, double, double> interactionTerm)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int k = 0; k < nz; k++)
                {
                    Matrix<double> matrix = Create.TridiagonalMatrixY(ny, ry);
                    Vector<double> vector = Create.RightHandSideY(solution, alpha, dt, ry, ny, i, k, dx, dy, dz, t);

                    // Добавляем член взаимодействия
                    for (int j = 0; j < ny; j++)
                    {
                        vector[j] += dt * interactionTerm(solution[i, j, k], otherSolution[i, j, k]);
                    }

                    Apply.BoundaryConditionsY(matrix, vector, ry, ny, dy, boundaryConditions, robinCoefficients);
                    Vector<double> newValues = matrix.Solve(vector);

                    for (int j = 0; j < nx; j++)
                    {
                        solution[i, j, k] = Math.Max(0, newValues[i]); // Гарантируем неотрицательность
                    }
                }
            }
        }

        private static void SolveZ(
        double[,,] solution,
        double alpha, double dt, double rz,
        int nx, int ny, int nz,
        double dx, double dy, double dz,
        Dictionary<string, (int type, double value)> boundaryConditions,
        Dictionary<string, (double alpha, double beta)> robinCoefficients,
        double t,
        double[,,] otherSolution,
        Func<double, double, double> interactionTerm)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    Matrix<double> matrix = Create.TridiagonalMatrixZ(nx, rz);
                    Vector<double> vector = Create.RightHandSideZ(solution, alpha, dt, rz, nx, i, j, dx, dy, dz, t);

                    // Добавляем член взаимодействия
                    for (int k = 0; k < nz; k++)
                    {
                        vector[k] += dt * interactionTerm(solution[i, j, k], otherSolution[i, j, k]);
                    }

                    Apply.BoundaryConditionsZ(matrix, vector, rz, nz, dz, boundaryConditions, robinCoefficients);
                    Vector<double> newValues = matrix.Solve(vector);

                    for (int k = 0; k < nz; k++)
                    {
                        solution[i, j, k] = Math.Max(0, newValues[i]); // Гарантируем неотрицательность
                    }
                }
            }
        }
    }
}

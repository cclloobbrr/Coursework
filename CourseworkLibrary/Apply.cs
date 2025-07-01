using MathNet.Numerics.LinearAlgebra;

namespace CourseworkLibrary
{
    public class Apply
    {
        public static void BoundaryConditionsX(Matrix<double> matrix, Vector<double> vector, double rx, int nx, int j, double dx, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients)
        {
            // Граница X0 
            if (boundaryConditions.ContainsKey("x0"))
            {
                var bc = boundaryConditions["x0"];
                switch (bc.type)
                {
                    case 0: // Дирихле
                        matrix.SetRow(0, Vector<double>.Build.Dense(nx));
                        matrix[0, 0] = 1;
                        vector[0] = bc.value;
                        break;
                    case 1: // Нейман
                        matrix[0, 0] = 1 + 2 * rx;
                        matrix[0, 1] = -2 * rx;
                        vector[0] = vector[0] + 2 * rx * dx * bc.value;
                        break;
                    case 2: // Робин
                        if (robinCoefficients == null || !robinCoefficients.ContainsKey("x0"))
                        {
                            throw new ArgumentException("Robin coefficients are missing for x0 boundary.");
                        }
                        var robin = robinCoefficients["x0"];
                        double alpha = robin.alpha;
                        double beta = robin.beta;

                        matrix[0, 0] = alpha * (1 + 2 * rx) + beta * (2 * rx / dx);
                        matrix[0, 1] = -alpha * (2 * rx);
                        vector[0] = alpha * vector[0] + 2 * beta * rx * bc.value;
                        break;
                }
            }
            else
            {
                Console.WriteLine("Boundary condition for x0 is missing.  Setting to zero gradient");
                matrix[0, 0] = 1 + 2 * rx;
                matrix[0, 1] = -2 * rx;
                vector[0] = vector[0];

            }

            // Граница X1
            if (boundaryConditions.ContainsKey("x1"))
            {
                var bc = boundaryConditions["x1"];
                switch (bc.type)
                {
                    case 0: // Дирихле
                        matrix.SetRow(nx - 1, Vector<double>.Build.Dense(nx));
                        matrix[nx - 1, nx - 1] = 1;
                        vector[nx - 1] = bc.value;
                        break;
                    case 1: // Нейман
                        matrix[nx - 1, nx - 1] = 1 + 2 * rx;
                        matrix[nx - 1, nx - 2] = -2 * rx;
                        vector[nx - 1] = vector[nx - 1] - 2 * rx * dx * bc.value;
                        break;
                    case 2: // Робин
                        if (robinCoefficients == null || !robinCoefficients.ContainsKey("x1"))
                        {
                            throw new ArgumentException("Robin coefficients are missing for x1 boundary.");
                        }
                        var robin = robinCoefficients["x1"];
                        double alpha = robin.alpha;
                        double beta = robin.beta;

                        matrix[nx - 1, nx - 1] = alpha * (1 + 2 * rx) + beta * (2 * rx / dx);
                        matrix[nx - 1, nx - 2] = -alpha * (2 * rx);

                        vector[nx - 1] = alpha * vector[nx - 1] + 2 * beta * rx * bc.value;
                        break;
                }
            }
            else
            {
                Console.WriteLine("Boundary condition for x1 is missing. Setting to zero gradient");
                matrix[nx - 1, nx - 1] = 1 + 2 * rx;
                matrix[nx - 1, nx - 2] = -2 * rx;
                vector[nx - 1] = vector[nx - 1];
            }
        }

        public static void BoundaryConditionsY(Matrix<double> matrix, Vector<double> vector, double ry, int ny, double dy, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients)
        {
            // Граница Y0
            if (boundaryConditions.ContainsKey("y0"))
            {
                var bc = boundaryConditions["y0"];
                switch (bc.type)
                {
                    case 0: // Дирихле
                        matrix.SetRow(0, Vector<double>.Build.Dense(ny));
                        matrix[0, 0] = 1;
                        vector[0] = bc.value;
                        break;
                    case 1: // Нейман
                        matrix[0, 0] = 1 + 2 * ry;
                        matrix[0, 1] = -2 * ry;
                        vector[0] = vector[0] + 2 * ry * dy * bc.value;
                        break;
                    case 2: // Робин
                        if (robinCoefficients == null || !robinCoefficients.ContainsKey("y0"))
                        {
                            throw new ArgumentException("Robin coefficients are missing for y0 boundary.");
                        }
                        var robin = robinCoefficients["y0"];
                        double alpha = robin.alpha;
                        double beta = robin.beta;
                        matrix[0, 0] = alpha * (1 + 2 * ry) + beta * (2 * ry / dy);
                        matrix[0, 1] = -alpha * (2 * ry);
                        vector[0] = alpha * vector[0] + 2 * beta * ry * bc.value;
                        break;
                }
            }
            else
            {
                Console.WriteLine("Boundary condition for y0 is missing.  Setting to zero gradient");
                matrix[0, 0] = 1 + 2 * ry;
                matrix[0, 1] = -2 * ry;
                vector[0] = vector[0];
            }

            // Граница Y1
            if (boundaryConditions.ContainsKey("y1"))
            {
                var bc = boundaryConditions["y1"];
                switch (bc.type)
                {
                    case 0: // Дирихле
                        matrix.SetRow(ny - 1, Vector<double>.Build.Dense(ny));
                        matrix[ny - 1, ny - 1] = 1;
                        vector[ny - 1] = bc.value;
                        break;
                    case 1: // Нейман
                        matrix[ny - 1, ny - 1] = 1 + 2 * ry;
                        matrix[ny - 1, ny - 2] = -2 * ry;
                        vector[ny - 1] = vector[ny - 1] - 2 * ry * dy * bc.value;
                        break;
                    case 2: // Робин
                        if (robinCoefficients == null || !robinCoefficients.ContainsKey("y1"))
                        {
                            throw new ArgumentException("Robin coefficients are missing for y1 boundary.");
                        }
                        var robin = robinCoefficients["y1"];
                        double alpha = robin.alpha;
                        double beta = robin.beta;
                        matrix[ny - 1, ny - 1] = alpha * (1 + 2 * ry) + beta * (2 * ry / dy);
                        matrix[ny - 1, ny - 2] = -alpha * (2 * ry);
                        vector[ny - 1] = alpha * vector[ny - 1] + 2 * beta * ry * bc.value;
                        break;
                }
            }
            else
            {
                Console.WriteLine("Boundary condition for y1 is missing. Setting to zero gradient");
                matrix[ny - 1, ny - 1] = 1 + 2 * ry;
                matrix[ny - 1, ny - 2] = -2 * ry;
                vector[ny - 1] = vector[ny - 1];

            }
        }

        public static void BoundaryConditionsZ(Matrix<double> matrix, Vector<double> vector, double rz, int nz, double dz, Dictionary<string, (int type, double value)> boundaryConditions, Dictionary<string, (double alpha, double beta)> robinCoefficients)
        {
            // Граница Z0 
            if (boundaryConditions.ContainsKey("z0"))
            {
                var bc = boundaryConditions["z0"];
                switch (bc.type)
                {
                    case 0: // Дирихле
                        matrix.SetRow(0, Vector<double>.Build.Dense(nz));
                        matrix[0, 0] = 1;
                        vector[0] = bc.value;
                        break;
                    case 1: // Нейман
                        matrix[0, 0] = 1 + 2 * rz;
                        matrix[0, 1] = -2 * rz;
                        vector[0] = vector[0] + 2 * rz * dz * bc.value;
                        break;
                    case 2: // Робин
                        if (robinCoefficients == null || !robinCoefficients.ContainsKey("z0"))
                        {
                            throw new ArgumentException("Robin coefficients are missing for z0 boundary.");
                        }
                        var robin = robinCoefficients["z0"];
                        double alpha = robin.alpha;
                        double beta = robin.beta;
                        matrix[0, 0] = alpha * (1 + 2 * rz) + beta * (2 * rz / dz);
                        matrix[0, 1] = -alpha * (2 * rz);
                        vector[0] = alpha * vector[0] + 2 * beta * rz * bc.value;
                        break;
                }
            }
            else
            {
                Console.WriteLine("Boundary condition for z0 is missing. Setting to zero gradient");
                matrix[0, 0] = 1 + 2 * rz;
                matrix[0, 1] = -2 * rz;
                vector[0] = vector[0];
            }

            // Граница Z1
            if (boundaryConditions.ContainsKey("z1"))
            {
                var bc = boundaryConditions["z1"];
                switch (bc.type)
                {
                    case 0: // Дирихле
                        matrix.SetRow(nz - 1, Vector<double>.Build.Dense(nz));
                        matrix[nz - 1, nz - 1] = 1;
                        vector[nz - 1] = bc.value;
                        break;
                    case 1: // Нейман
                        matrix[nz - 1, nz - 1] = 1 + 2 * rz;
                        matrix[nz - 1, nz - 2] = -2 * rz;
                        vector[nz - 1] = vector[nz - 1] - 2 * rz * dz * bc.value;
                        break;
                    case 2: // Робин
                        if (robinCoefficients == null || !robinCoefficients.ContainsKey("z1"))
                        {
                            throw new ArgumentException("Robin coefficients are missing for z1 boundary.");
                        }
                        var robin = robinCoefficients["z1"];
                        double alpha = robin.alpha;
                        double beta = robin.beta;
                        matrix[nz - 1, nz - 1] = alpha * (1 + 2 * rz) + beta * (2 * rz / dz);
                        matrix[nz - 1, nz - 2] = -alpha * (2 * rz);
                        vector[nz - 1] = alpha * vector[nz - 1] + 2 * beta * rz * bc.value;
                        break;
                }
            }
            else
            {
                Console.WriteLine("Boundary condition for z1 is missing. Setting to zero gradient");
                matrix[nz - 1, nz - 1] = 1 + 2 * rz;
                matrix[nz - 1, nz - 2] = -2 * rz;
                vector[nz - 1] = vector[nz - 1];
            }
        }
    }
}
